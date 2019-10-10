% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.18s
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
	% StartTime: 2019-10-10 11:20:31
	% EndTime: 2019-10-10 11:20:31
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (179->35), mult. (217->48), div. (0->0), fcn. (148->6), ass. (0->34)
	t183 = qJ(2) + qJ(3);
	t181 = cos(t183);
	t211 = r_i_i_C(3) + qJ(4);
	t195 = t211 * t181;
	t180 = sin(t183);
	t178 = t180 * qJD(4);
	t182 = qJD(2) + qJD(3);
	t214 = pkin(3) + r_i_i_C(1);
	t203 = t214 * t180;
	t184 = sin(qJ(2));
	t210 = pkin(2) * qJD(2);
	t204 = t184 * t210;
	t217 = (-t203 + t195) * t182 + (r_i_i_C(2) + pkin(8) + pkin(7)) * qJD(1) + t178 - t204;
	t213 = pkin(2) * t184;
	t209 = t181 * t182;
	t187 = cos(qJ(1));
	t208 = t182 * t187;
	t185 = sin(qJ(1));
	t207 = qJD(1) * t185;
	t206 = qJD(1) * t187;
	t205 = qJD(4) * t181;
	t202 = t214 * t187;
	t201 = t185 * t209;
	t200 = t185 * t205 + t206 * t195;
	t197 = t180 * t207;
	t199 = t187 * t205 + t214 * t197;
	t196 = t211 * t180;
	t194 = t211 * t185;
	t192 = -t214 * t181 - t196;
	t191 = -t182 * t203 + t211 * t209 + t178;
	t186 = cos(qJ(2));
	t190 = qJD(1) * (-t186 * pkin(2) - pkin(1) + t192);
	t189 = t192 * t182 - t186 * t210;
	t1 = [-t217 * t185 + t187 * t190, (-t195 + t213) * t207 + t189 * t187 + t199, -t196 * t208 + (-qJD(1) * t194 - t182 * t202) * t181 + t199, t181 * t208 - t197, 0, 0; t185 * t190 + t217 * t187, (-t203 - t213) * t206 + t189 * t185 + t200, -t214 * t201 + (-qJD(1) * t202 - t182 * t194) * t180 + t200, t180 * t206 + t201, 0, 0; 0, t191 - t204, t191, t182 * t180, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (221->46), mult. (265->57), div. (0->0), fcn. (180->6), ass. (0->36)
	t68 = qJ(2) + qJ(3);
	t66 = cos(t68);
	t97 = r_i_i_C(1) + qJ(4);
	t109 = t97 * t66;
	t65 = sin(t68);
	t63 = t65 * qJD(4);
	t67 = qJD(2) + qJD(3);
	t69 = sin(qJ(2));
	t96 = pkin(2) * qJD(2);
	t88 = t69 * t96;
	t110 = pkin(3) + pkin(4);
	t92 = -r_i_i_C(2) + t110;
	t111 = (-t92 * t65 + t109) * t67 - (r_i_i_C(3) + qJ(5) - pkin(8) - pkin(7)) * qJD(1) + t63 - t88;
	t70 = sin(qJ(1));
	t100 = t66 * t70;
	t72 = cos(qJ(1));
	t94 = qJD(1) * t72;
	t108 = t67 * t100 + t65 * t94;
	t102 = pkin(2) * t69;
	t99 = t67 * t65;
	t98 = t67 * t72;
	t95 = qJD(1) * t70;
	t93 = qJD(4) * t66;
	t90 = t66 * t98;
	t87 = t110 * t67;
	t86 = t110 * t72;
	t84 = t65 * t95;
	t82 = t97 * t65;
	t80 = t97 * t70;
	t79 = r_i_i_C(2) * t90 + t110 * t84 + t72 * t93;
	t78 = t108 * r_i_i_C(2) + t109 * t94 + t70 * t93;
	t76 = r_i_i_C(2) * t99 + t109 * t67 - t65 * t87 + t63;
	t71 = cos(qJ(2));
	t75 = -qJD(5) + (-t71 * pkin(2) - t92 * t66 - pkin(1) - t82) * qJD(1);
	t74 = -t71 * t96 + (-t110 * t66 - t82) * t67;
	t1 = [-t111 * t70 + t75 * t72, (-r_i_i_C(2) * t65 + t102 - t109) * t95 + t74 * t72 + t79, (-r_i_i_C(2) * t95 - t97 * t98) * t65 + (-qJD(1) * t80 - t67 * t86) * t66 + t79, -t84 + t90, -t94, 0; t111 * t72 + t75 * t70, (-t110 * t65 - t102) * t94 + t74 * t70 + t78, -t87 * t100 + (-qJD(1) * t86 - t67 * t80) * t65 + t78, t108, -t95, 0; 0, t76 - t88, t76, t99, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:31
	% EndTime: 2019-10-10 11:20:32
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (410->65), mult. (545->89), div. (0->0), fcn. (404->8), ass. (0->55)
	t264 = qJ(2) + qJ(3);
	t262 = cos(t264);
	t268 = cos(qJ(6));
	t311 = pkin(5) + qJ(4);
	t319 = r_i_i_C(1) * t268 + t311;
	t320 = t262 * t319;
	t261 = sin(t264);
	t259 = t261 * qJD(4);
	t263 = qJD(2) + qJD(3);
	t297 = pkin(3) + pkin(4) + pkin(9) + r_i_i_C(3);
	t286 = t297 * t261;
	t266 = sin(qJ(2));
	t310 = pkin(2) * qJD(2);
	t298 = t266 * t310;
	t318 = (-t311 * t262 + t286) * t263 + (qJ(5) - pkin(8) - pkin(7)) * qJD(1) - t259 + t298;
	t265 = sin(qJ(6));
	t308 = t263 * t261;
	t296 = t265 * t308;
	t317 = r_i_i_C(2) * t296 + qJD(4) * t262;
	t313 = pkin(2) * t266;
	t309 = t262 * t263;
	t307 = t263 * t268;
	t270 = cos(qJ(1));
	t306 = t263 * t270;
	t305 = t268 * t270;
	t267 = sin(qJ(1));
	t303 = qJD(1) * t267;
	t302 = qJD(1) * t270;
	t300 = qJD(6) * t262;
	t299 = r_i_i_C(2) * t262 * t265;
	t295 = t262 * t307;
	t294 = t267 * t309;
	t293 = t262 * t306;
	t291 = t261 * t303;
	t289 = qJD(6) * t261 + qJD(1);
	t288 = qJD(1) * t261 + qJD(6);
	t285 = t297 * t262;
	t284 = t317 * t267 + t302 * t320;
	t283 = t289 * t265;
	t282 = t319 * t261;
	t280 = (-r_i_i_C(1) * t265 - r_i_i_C(2) * t268) * qJD(6);
	t279 = t317 * t270 + t297 * t291 + t299 * t303;
	t278 = -t286 - t299;
	t277 = t288 * t267 - t293;
	t276 = -t297 * t263 + t280;
	t269 = cos(qJ(2));
	t275 = -qJD(5) + (-t269 * pkin(2) - t311 * t261 - pkin(1) - t285) * qJD(1);
	t274 = t262 * t280 + (-t282 - t285) * t263;
	t273 = r_i_i_C(1) * t295 + t276 * t261 - t263 * t299 + t311 * t309 + t259;
	t272 = -t269 * t310 + t274;
	t243 = -t288 * t305 + (t283 - t295) * t267;
	t242 = t289 * t268 * t267 + (t288 * t270 + t294) * t265;
	t241 = t277 * t268 + t270 * t283;
	t240 = t277 * t265 - t289 * t305;
	t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) + t318 * t267 + t275 * t270, (t313 - t320) * t303 + t272 * t270 + t279, -t282 * t306 + (t276 * t270 - t303 * t319) * t262 + t279, -t291 + t293, -t302, t240 * r_i_i_C(1) + t241 * r_i_i_C(2); -t241 * r_i_i_C(1) + t240 * r_i_i_C(2) + t275 * t267 - t318 * t270, (t278 - t313) * t302 + t272 * t267 + t284, t274 * t267 + t278 * t302 + t284, t261 * t302 + t294, -t303, -t242 * r_i_i_C(1) + t243 * r_i_i_C(2); 0, t273 - t298, t273, t308, 0, (-t261 * t307 - t265 * t300) * r_i_i_C(2) + (t268 * t300 - t296) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end