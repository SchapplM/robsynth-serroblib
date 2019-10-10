% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
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
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.10s
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(8) + pkin(7);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (149->33), mult. (146->42), div. (0->0), fcn. (95->8), ass. (0->32)
	t52 = sin(qJ(2));
	t50 = qJD(2) + qJD(3);
	t51 = qJ(2) + qJ(3);
	t46 = pkin(10) + t51;
	t44 = sin(t46);
	t45 = cos(t46);
	t59 = r_i_i_C(1) * t44 + r_i_i_C(2) * t45;
	t47 = sin(t51);
	t75 = pkin(3) * t47;
	t78 = t59 + t75;
	t56 = t78 * t50;
	t67 = pkin(2) * qJD(2);
	t76 = t52 * t67 + t56;
	t48 = cos(t51);
	t72 = r_i_i_C(1) * t45;
	t77 = -pkin(3) * t48 - t72;
	t71 = r_i_i_C(2) * t44;
	t69 = r_i_i_C(3) + qJ(4) + pkin(8) + pkin(7);
	t68 = t48 * t50;
	t53 = sin(qJ(1));
	t66 = qJD(1) * t53;
	t55 = cos(qJ(1));
	t65 = qJD(1) * t55;
	t64 = t50 * t72;
	t63 = t50 * t71;
	t62 = t55 * t63 + t59 * t66;
	t54 = cos(qJ(2));
	t60 = -pkin(3) * t68 - t54 * t67 - t64;
	t58 = -t54 * pkin(2) - pkin(1) + t71 + t77;
	t43 = -t52 * pkin(2) - t75;
	t38 = t53 * t63;
	t1 = [t55 * qJD(4) + t76 * t53 + (-t69 * t53 + t58 * t55) * qJD(1), -t43 * t66 + t60 * t55 + t62, -t55 * t64 + (t47 * t66 - t55 * t68) * pkin(3) + t62, t65, 0, 0; t53 * qJD(4) - t76 * t55 + (t58 * t53 + t69 * t55) * qJD(1), t38 + t60 * t53 + (t43 - t59) * t65, t50 * t53 * t77 - t65 * t78 + t38, t66, 0, 0; 0, -t76, -t56, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:49
	% EndTime: 2019-10-10 11:18:49
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (288->48), mult. (249->53), div. (0->0), fcn. (170->8), ass. (0->42)
	t201 = qJ(2) + qJ(3);
	t196 = pkin(10) + t201;
	t195 = cos(t196);
	t230 = r_i_i_C(3) + qJ(5);
	t215 = t230 * t195;
	t194 = sin(t196);
	t193 = t194 * qJD(5);
	t197 = sin(t201);
	t200 = qJD(2) + qJD(3);
	t202 = sin(qJ(2));
	t229 = pkin(2) * qJD(2);
	t223 = t202 * t229;
	t233 = pkin(3) * t200;
	t236 = pkin(4) - r_i_i_C(2);
	t241 = (-t236 * t194 + t215) * t200 + (r_i_i_C(1) + qJ(4) + pkin(8) + pkin(7)) * qJD(1) - t197 * t233 + t193 - t223;
	t203 = sin(qJ(1));
	t205 = cos(qJ(1));
	t225 = qJD(1) * t205;
	t228 = t195 * t200;
	t240 = t194 * t225 + t203 * t228;
	t235 = pkin(3) * t197;
	t198 = cos(t201);
	t234 = pkin(3) * t198;
	t232 = pkin(4) * t194;
	t227 = t200 * t194;
	t226 = qJD(1) * t203;
	t224 = qJD(5) * t195;
	t221 = t205 * t228;
	t218 = t194 * t226;
	t220 = pkin(4) * t218 + r_i_i_C(2) * t221 + t205 * t224;
	t216 = t230 * t194;
	t213 = -t232 - t235;
	t212 = t240 * r_i_i_C(2) + t203 * t224 + t225 * t215;
	t211 = -pkin(4) * t195 - t216;
	t210 = -r_i_i_C(2) * t194 - t215;
	t204 = cos(qJ(2));
	t209 = -t198 * t233 + t211 * t200 - t204 * t229;
	t208 = t200 * (t211 - t234);
	t207 = r_i_i_C(2) * t227 + t213 * t200 + t230 * t228 + t193;
	t206 = qJD(4) + (-t204 * pkin(2) - t236 * t195 - pkin(1) - t216 - t234) * qJD(1);
	t190 = -t202 * pkin(2) - t235;
	t1 = [-t241 * t203 + t206 * t205, t209 * t205 + (-t190 + t210) * t226 + t220, t205 * t208 + (t210 + t235) * t226 + t220, t225, -t218 + t221, 0; t206 * t203 + t241 * t205, (t190 - t232) * t225 + t209 * t203 + t212, t203 * t208 + t213 * t225 + t212, t226, t240, 0; 0, t207 - t223, t207, 0, t227, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:49
	% EndTime: 2019-10-10 11:18:50
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (503->69), mult. (493->90), div. (0->0), fcn. (370->10), ass. (0->60)
	t272 = qJ(2) + qJ(3);
	t267 = pkin(10) + t272;
	t266 = cos(t267);
	t273 = sin(qJ(6));
	t276 = cos(qJ(6));
	t327 = r_i_i_C(1) * t273 + r_i_i_C(2) * t276 + qJ(5);
	t285 = t327 * t266;
	t265 = sin(t267);
	t274 = sin(qJ(2));
	t316 = pkin(2) * qJD(2);
	t305 = t274 * t316;
	t268 = sin(t272);
	t271 = qJD(2) + qJD(3);
	t315 = t266 * t271;
	t320 = pkin(3) * t271;
	t324 = -qJ(5) * t315 + t268 * t320;
	t307 = pkin(4) + pkin(9) + r_i_i_C(3);
	t328 = -t307 * t271 + qJD(5);
	t330 = t328 * t265 + (pkin(5) + qJ(4) + pkin(8) + pkin(7)) * qJD(1) - t305 - t324;
	t308 = qJD(6) * t276;
	t298 = t266 * t308;
	t329 = r_i_i_C(1) * t298 + qJD(5) * t266;
	t292 = qJD(6) * t265 + qJD(1);
	t326 = t276 * t292;
	t325 = t292 * t273;
	t322 = pkin(3) * t268;
	t269 = cos(t272);
	t321 = pkin(3) * t269;
	t314 = t271 * t273;
	t313 = t271 * t276;
	t275 = sin(qJ(1));
	t312 = qJD(1) * t275;
	t278 = cos(qJ(1));
	t311 = qJD(1) * t278;
	t309 = qJD(6) * t273;
	t304 = t266 * t313;
	t303 = t275 * t315;
	t302 = t278 * t315;
	t300 = t265 * t312;
	t299 = t266 * t309;
	t297 = t307 * t265;
	t296 = t307 * t266;
	t293 = r_i_i_C(2) * t299;
	t291 = -qJD(1) * t265 - qJD(6);
	t289 = t329 * t275 + t311 * t285;
	t288 = t329 * t278 + t307 * t300;
	t287 = t291 * t278;
	t284 = t291 * t275 + t302;
	t277 = cos(qJ(2));
	t283 = qJD(4) + (-t277 * pkin(2) - qJ(5) * t265 - pkin(1) - t296 - t321) * qJD(1);
	t282 = -t265 * t327 - t296;
	t281 = t266 * r_i_i_C(1) * t314 + r_i_i_C(2) * t304 - t324 + (r_i_i_C(1) * t308 - r_i_i_C(2) * t309 + t328) * t265;
	t280 = -t269 * t320 + t282 * t271 - t277 * t316 - t293;
	t279 = -t293 + (t282 - t321) * t271;
	t261 = -t274 * pkin(2) - t322;
	t245 = t284 * t273 + t278 * t326;
	t244 = t284 * t276 - t278 * t325;
	t243 = -t275 * t326 + (t287 - t303) * t273;
	t242 = t276 * t287 + (-t304 + t325) * t275;
	t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) - t330 * t275 + t283 * t278, (-t261 - t285) * t312 + t280 * t278 + t288, (-t285 + t322) * t312 + t279 * t278 + t288, t311, -t300 + t302, t244 * r_i_i_C(1) - t245 * r_i_i_C(2); t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t283 * t275 + t330 * t278, (t261 - t297) * t311 + t280 * t275 + t289, (-t297 - t322) * t311 + t279 * t275 + t289, t312, t265 * t311 + t303, -t242 * r_i_i_C(1) + t243 * r_i_i_C(2); 0, t281 - t305, t281, 0, t271 * t265, (-t265 * t314 + t298) * r_i_i_C(2) + (t265 * t313 + t299) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end