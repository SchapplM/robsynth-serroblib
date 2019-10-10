% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:06
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t20 * t22) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t18 * t22 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(10);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(7);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0, 0; 0, -t33, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (131->30), mult. (132->40), div. (0->0), fcn. (86->8), ass. (0->27)
	t49 = qJ(2) + pkin(10);
	t41 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t49);
	t48 = qJD(2) + qJD(4);
	t46 = qJ(4) + t49;
	t43 = cos(t46);
	t67 = r_i_i_C(2) * t43;
	t42 = sin(t46);
	t69 = r_i_i_C(1) * t42;
	t56 = t67 + t69;
	t54 = t56 * t48;
	t70 = t41 * qJD(2) - t54;
	t68 = r_i_i_C(2) * t42;
	t66 = r_i_i_C(3) + pkin(8) + qJ(3) + pkin(7);
	t65 = t43 * t48;
	t51 = sin(qJ(1));
	t64 = qJD(1) * t51;
	t53 = cos(qJ(1));
	t63 = qJD(1) * t53;
	t62 = r_i_i_C(1) * t65;
	t61 = t48 * t68;
	t59 = qJD(1) * t67;
	t60 = t51 * t59 + t53 * t61 + t64 * t69;
	t57 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t49);
	t58 = t57 * qJD(2) - t62;
	t55 = -r_i_i_C(1) * t43 - pkin(1) + t57 + t68;
	t36 = t51 * t61;
	t1 = [t53 * qJD(3) - t70 * t51 + (-t66 * t51 + t55 * t53) * qJD(1), -t41 * t64 + t58 * t53 + t60, t63, -t53 * t62 + t60, 0, 0; t51 * qJD(3) + t70 * t53 + (t55 * t51 + t66 * t53) * qJD(1), t36 + t58 * t51 + (t41 - t56) * t63, t64, -t53 * t59 + t36 + (-t42 * t63 - t51 * t65) * r_i_i_C(1), 0, 0; 0, t70, 0, -t54, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:10
	% EndTime: 2019-10-10 10:06:10
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (270->43), mult. (235->50), div. (0->0), fcn. (161->8), ass. (0->34)
	t199 = qJ(2) + pkin(10);
	t196 = qJ(4) + t199;
	t193 = cos(t196);
	t226 = r_i_i_C(3) + qJ(5);
	t212 = t226 * t193;
	t186 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t199);
	t176 = t186 * qJD(2);
	t192 = sin(t196);
	t191 = t192 * qJD(5);
	t198 = qJD(2) + qJD(4);
	t228 = pkin(4) - r_i_i_C(2);
	t233 = (-t228 * t192 + t212) * t198 + (r_i_i_C(1) + pkin(8) + qJ(3) + pkin(7)) * qJD(1) + t176 + t191;
	t201 = sin(qJ(1));
	t223 = t198 * t201;
	t219 = t193 * t223;
	t203 = cos(qJ(1));
	t221 = qJD(1) * t203;
	t232 = t192 * t221 + t219;
	t225 = t193 * t198;
	t224 = t198 * t192;
	t222 = qJD(1) * t201;
	t220 = qJD(5) * t193;
	t218 = t203 * t225;
	t215 = t192 * t222;
	t217 = pkin(4) * t215 + r_i_i_C(2) * t218 + t203 * t220;
	t213 = t226 * t192;
	t210 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t199);
	t209 = t232 * r_i_i_C(2) + t201 * t220 + t221 * t212;
	t208 = -r_i_i_C(2) * t192 - t212;
	t207 = -t228 * t224 + t226 * t225 + t191;
	t206 = (-pkin(4) * t193 - t213) * t198;
	t205 = t210 * qJD(2) + t206;
	t204 = qJD(3) + (-t228 * t193 - pkin(1) + t210 - t213) * qJD(1);
	t1 = [-t233 * t201 + t204 * t203, t205 * t203 + (-t186 + t208) * t222 + t217, t221, t203 * t206 + t208 * t222 + t217, -t215 + t218, 0; t204 * t201 + t233 * t203, (-pkin(4) * t192 + t186) * t221 + t205 * t201 + t209, t222, -pkin(4) * t219 + (-pkin(4) * t221 - t226 * t223) * t192 + t209, t232, 0; 0, t176 + t207, 0, t207, t224, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:10
	% EndTime: 2019-10-10 10:06:10
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (485->65), mult. (479->88), div. (0->0), fcn. (361->10), ass. (0->54)
	t268 = qJ(2) + pkin(10);
	t265 = qJ(4) + t268;
	t262 = cos(t265);
	t269 = sin(qJ(6));
	t272 = cos(qJ(6));
	t317 = r_i_i_C(1) * t269 + r_i_i_C(2) * t272 + qJ(5);
	t321 = t262 * t317;
	t257 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t268);
	t246 = t257 * qJD(2);
	t267 = qJD(2) + qJD(4);
	t310 = t262 * t267;
	t256 = qJ(5) * t310;
	t261 = sin(t265);
	t301 = pkin(4) + pkin(9) + r_i_i_C(3);
	t291 = t301 * t267;
	t320 = (qJD(5) - t291) * t261 + (pkin(5) + pkin(8) + qJ(3) + pkin(7)) * qJD(1) + t246 + t256;
	t302 = qJD(6) * t272;
	t294 = t262 * t302;
	t319 = r_i_i_C(1) * t294 + qJD(5) * t262;
	t289 = qJD(6) * t261 + qJD(1);
	t316 = t272 * t289;
	t315 = t289 * t269;
	t309 = t267 * t269;
	t308 = t267 * t272;
	t274 = cos(qJ(1));
	t307 = t267 * t274;
	t271 = sin(qJ(1));
	t306 = qJD(1) * t271;
	t305 = qJD(1) * t274;
	t303 = qJD(6) * t269;
	t300 = t262 * t308;
	t299 = t271 * t310;
	t298 = t262 * t307;
	t296 = t261 * t306;
	t295 = t262 * t303;
	t293 = t301 * t261;
	t292 = t301 * t262;
	t288 = -qJD(1) * t261 - qJD(6);
	t287 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t268);
	t285 = t319 * t271 + t305 * t321;
	t284 = t319 * t274 + t301 * t296;
	t283 = t288 * t274;
	t281 = t317 * t261;
	t280 = -r_i_i_C(2) * t303 - t291;
	t279 = t288 * t271 + t298;
	t278 = qJD(3) + (-qJ(5) * t261 - pkin(1) + t287 - t292) * qJD(1);
	t277 = t262 * r_i_i_C(1) * t309 + r_i_i_C(2) * t300 + t256 + (r_i_i_C(1) * t302 + qJD(5) + t280) * t261;
	t276 = -r_i_i_C(2) * t295 + (-t292 - t281) * t267;
	t275 = t287 * qJD(2) + t276;
	t241 = t279 * t269 + t274 * t316;
	t240 = t279 * t272 - t274 * t315;
	t239 = -t271 * t316 + (t283 - t299) * t269;
	t238 = t272 * t283 + (-t300 + t315) * t271;
	t1 = [t239 * r_i_i_C(1) + t238 * r_i_i_C(2) - t320 * t271 + t278 * t274, (-t257 - t321) * t306 + t275 * t274 + t284, t305, -t281 * t307 + (t280 * t274 - t306 * t317) * t262 + t284, -t296 + t298, t240 * r_i_i_C(1) - t241 * r_i_i_C(2); t241 * r_i_i_C(1) + t240 * r_i_i_C(2) + t278 * t271 + t320 * t274, (t257 - t293) * t305 + t275 * t271 + t285, t306, t276 * t271 - t293 * t305 + t285, t261 * t305 + t299, -t238 * r_i_i_C(1) + t239 * r_i_i_C(2); 0, t246 + t277, 0, t277, t267 * t261, (-t261 * t309 + t294) * r_i_i_C(2) + (t261 * t308 + t295) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end