% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(10)) + r_i_i_C(2) * sin(pkin(10)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(10) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(10)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (124->30), mult. (120->40), div. (0->0), fcn. (79->7), ass. (0->28)
	t44 = qJD(3) + qJD(4);
	t43 = pkin(10) + qJ(3);
	t41 = qJ(4) + t43;
	t38 = cos(t41);
	t61 = r_i_i_C(2) * t38;
	t37 = sin(t41);
	t63 = r_i_i_C(1) * t37;
	t51 = t61 + t63;
	t49 = t51 * t44;
	t39 = sin(t43);
	t64 = pkin(3) * t39;
	t65 = qJD(3) * t64 + t49;
	t62 = r_i_i_C(2) * t37;
	t60 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(2);
	t59 = t38 * t44;
	t45 = sin(qJ(1));
	t58 = qJD(1) * t45;
	t46 = cos(qJ(1));
	t57 = qJD(1) * t46;
	t40 = cos(t43);
	t56 = qJD(3) * t40;
	t55 = r_i_i_C(1) * t59;
	t54 = t44 * t62;
	t52 = qJD(1) * t61;
	t50 = -r_i_i_C(1) * t38 - pkin(3) * t40 - cos(pkin(10)) * pkin(2) - pkin(1) + t62;
	t48 = t45 * t52 + t58 * t63 + (t54 - t55) * t46;
	t32 = t45 * t54;
	t1 = [t46 * qJD(2) + t65 * t45 + (-t60 * t45 + t50 * t46) * qJD(1), t57, (t39 * t58 - t46 * t56) * pkin(3) + t48, t48, 0, 0; t45 * qJD(2) - t65 * t46 + (t50 * t45 + t60 * t46) * qJD(1), t58, t32 + (-pkin(3) * t56 - t55) * t45 + (-t51 - t64) * t57, -t46 * t52 + t32 + (-t37 * t57 - t45 * t59) * r_i_i_C(1), 0, 0; 0, 0, -t65, -t49, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:30
	% EndTime: 2019-10-10 01:30:31
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (263->42), mult. (223->50), div. (0->0), fcn. (154->7), ass. (0->36)
	t193 = pkin(10) + qJ(3);
	t191 = qJ(4) + t193;
	t188 = cos(t191);
	t220 = r_i_i_C(3) + qJ(5);
	t204 = t220 * t188;
	t187 = sin(t191);
	t186 = t187 * qJD(5);
	t194 = qJD(3) + qJD(4);
	t189 = sin(t193);
	t219 = pkin(3) * qJD(3);
	t212 = t189 * t219;
	t223 = pkin(4) - r_i_i_C(2);
	t228 = (-t223 * t187 + t204) * t194 + (r_i_i_C(1) + pkin(8) + pkin(7) + qJ(2)) * qJD(1) + t186 - t212;
	t195 = sin(qJ(1));
	t216 = t194 * t195;
	t211 = t188 * t216;
	t196 = cos(qJ(1));
	t214 = qJD(1) * t196;
	t227 = t187 * t214 + t211;
	t222 = pkin(3) * t189;
	t218 = t188 * t194;
	t217 = t194 * t187;
	t215 = qJD(1) * t195;
	t213 = qJD(5) * t188;
	t210 = t196 * t218;
	t207 = t187 * t215;
	t209 = pkin(4) * t207 + r_i_i_C(2) * t210 + t196 * t213;
	t205 = t220 * t187;
	t203 = t227 * r_i_i_C(2) + t195 * t213 + t214 * t204;
	t202 = -r_i_i_C(2) * t187 - t204;
	t200 = -t223 * t217 + t220 * t218 + t186;
	t199 = (-pkin(4) * t188 - t205) * t194;
	t190 = cos(t193);
	t198 = -t190 * t219 + t199;
	t197 = qJD(2) + (-t223 * t188 - pkin(3) * t190 - cos(pkin(10)) * pkin(2) - pkin(1) - t205) * qJD(1);
	t1 = [-t228 * t195 + t197 * t196, t214, t198 * t196 + (t202 + t222) * t215 + t209, t196 * t199 + t202 * t215 + t209, -t207 + t210, 0; t197 * t195 + t228 * t196, t215, (-pkin(4) * t187 - t222) * t214 + t198 * t195 + t203, -pkin(4) * t211 + (-pkin(4) * t214 - t220 * t216) * t187 + t203, t227, 0; 0, 0, t200 - t212, t200, t217, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:31
	% EndTime: 2019-10-10 01:30:31
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (478->64), mult. (467->88), div. (0->0), fcn. (354->9), ass. (0->56)
	t262 = pkin(10) + qJ(3);
	t260 = qJ(4) + t262;
	t257 = cos(t260);
	t264 = sin(qJ(6));
	t266 = cos(qJ(6));
	t312 = r_i_i_C(1) * t264 + r_i_i_C(2) * t266 + qJ(5);
	t316 = t257 * t312;
	t263 = qJD(3) + qJD(4);
	t303 = t257 * t263;
	t252 = qJ(5) * t303;
	t256 = sin(t260);
	t294 = pkin(4) + pkin(9) + r_i_i_C(3);
	t283 = t294 * t263;
	t258 = sin(t262);
	t304 = pkin(3) * qJD(3);
	t293 = t258 * t304;
	t315 = (qJD(5) - t283) * t256 + (pkin(5) + pkin(8) + pkin(7) + qJ(2)) * qJD(1) + t252 - t293;
	t295 = qJD(6) * t266;
	t286 = t257 * t295;
	t314 = r_i_i_C(1) * t286 + qJD(5) * t257;
	t281 = qJD(6) * t256 + qJD(1);
	t311 = t266 * t281;
	t310 = t281 * t264;
	t308 = pkin(3) * t258;
	t302 = t263 * t264;
	t301 = t263 * t266;
	t267 = cos(qJ(1));
	t300 = t263 * t267;
	t265 = sin(qJ(1));
	t299 = qJD(1) * t265;
	t298 = qJD(1) * t267;
	t296 = qJD(6) * t264;
	t292 = t257 * t301;
	t291 = t265 * t303;
	t290 = t257 * t300;
	t288 = t256 * t299;
	t287 = t257 * t296;
	t285 = t294 * t256;
	t284 = t294 * t257;
	t280 = -qJD(1) * t256 - qJD(6);
	t279 = t314 * t265 + t298 * t316;
	t278 = t314 * t267 + t294 * t288;
	t277 = t280 * t267;
	t274 = t312 * t256;
	t273 = -r_i_i_C(2) * t296 - t283;
	t272 = t280 * t265 + t290;
	t259 = cos(t262);
	t271 = qJD(2) + (-qJ(5) * t256 - pkin(3) * t259 - cos(pkin(10)) * pkin(2) - pkin(1) - t284) * qJD(1);
	t270 = t257 * r_i_i_C(1) * t302 + r_i_i_C(2) * t292 + t252 + (r_i_i_C(1) * t295 + qJD(5) + t273) * t256;
	t269 = -r_i_i_C(2) * t287 + (-t284 - t274) * t263;
	t268 = -t259 * t304 + t269;
	t239 = t272 * t264 + t267 * t311;
	t238 = t272 * t266 - t267 * t310;
	t237 = -t265 * t311 + (t277 - t291) * t264;
	t236 = t266 * t277 + (-t292 + t310) * t265;
	t1 = [t237 * r_i_i_C(1) + t236 * r_i_i_C(2) - t315 * t265 + t271 * t267, t298, (t308 - t316) * t299 + t268 * t267 + t278, -t274 * t300 + (t273 * t267 - t299 * t312) * t257 + t278, -t288 + t290, t238 * r_i_i_C(1) - t239 * r_i_i_C(2); t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t271 * t265 + t315 * t267, t299, (-t285 - t308) * t298 + t268 * t265 + t279, t269 * t265 - t285 * t298 + t279, t256 * t298 + t291, -t236 * r_i_i_C(1) + t237 * r_i_i_C(2); 0, 0, t270 - t293, t270, t263 * t256, (-t256 * t302 + t286) * r_i_i_C(2) + (t256 * t301 + t287) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end