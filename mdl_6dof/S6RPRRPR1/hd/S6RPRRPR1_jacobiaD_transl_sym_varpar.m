% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
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
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(10);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(7) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (119->29), mult. (118->40), div. (0->0), fcn. (75->8), ass. (0->26)
	t44 = qJD(3) + qJD(4);
	t46 = qJ(3) + qJ(4);
	t42 = sin(t46);
	t43 = cos(t46);
	t63 = r_i_i_C(2) * t43;
	t54 = r_i_i_C(1) * t42 + t63;
	t52 = t54 * t44;
	t47 = sin(qJ(3));
	t65 = pkin(3) * t47;
	t66 = qJD(3) * t65 + t52;
	t64 = r_i_i_C(2) * t42;
	t62 = r_i_i_C(3) + pkin(8) + pkin(7);
	t61 = t43 * t44;
	t60 = qJD(1) * t42;
	t48 = cos(qJ(3));
	t59 = qJD(3) * t48;
	t58 = r_i_i_C(1) * t61;
	t57 = t44 * t64;
	t56 = qJD(1) * t63;
	t53 = -t48 * pkin(3) - r_i_i_C(1) * t43 - pkin(2) + t64;
	t45 = qJ(1) + pkin(10);
	t40 = sin(t45);
	t41 = cos(t45);
	t51 = (t57 - t58) * t41 + (t60 * r_i_i_C(1) + t56) * t40;
	t35 = t40 * t57;
	t1 = [t66 * t40 + (-cos(qJ(1)) * pkin(1) - t62 * t40 + t53 * t41) * qJD(1), 0, (qJD(1) * t40 * t47 - t41 * t59) * pkin(3) + t51, t51, 0, 0; -t66 * t41 + (-sin(qJ(1)) * pkin(1) + t62 * t41 + t53 * t40) * qJD(1), 0, t35 + (-pkin(3) * t59 - t58) * t40 + (-t54 - t65) * t41 * qJD(1), -t41 * t56 + t35 + (-t40 * t61 - t41 * t60) * r_i_i_C(1), 0, 0; 0, 0, -t66, -t52, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:32
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (195->36), mult. (150->44), div. (0->0), fcn. (97->10), ass. (0->33)
	t59 = sin(qJ(3));
	t56 = qJD(3) + qJD(4);
	t58 = qJ(3) + qJ(4);
	t52 = pkin(11) + t58;
	t48 = sin(t52);
	t49 = cos(t52);
	t64 = r_i_i_C(1) * t48 + r_i_i_C(2) * t49;
	t53 = sin(t58);
	t80 = pkin(4) * t53;
	t83 = t64 + t80;
	t61 = t83 * t56;
	t72 = pkin(3) * qJD(3);
	t81 = t59 * t72 + t61;
	t54 = cos(t58);
	t77 = r_i_i_C(1) * t49;
	t82 = -pkin(4) * t54 - t77;
	t76 = r_i_i_C(2) * t48;
	t74 = r_i_i_C(3) + qJ(5) + pkin(8) + pkin(7);
	t73 = t54 * t56;
	t57 = qJ(1) + pkin(10);
	t50 = sin(t57);
	t71 = qJD(1) * t50;
	t51 = cos(t57);
	t70 = qJD(1) * t51;
	t69 = t56 * t77;
	t68 = t56 * t76;
	t67 = t51 * t68 + t64 * t71;
	t60 = cos(qJ(3));
	t65 = -pkin(4) * t73 - t60 * t72 - t69;
	t63 = -pkin(3) * t60 - pkin(2) + t76 + t82;
	t47 = -pkin(3) * t59 - t80;
	t40 = t50 * t68;
	t1 = [t51 * qJD(5) + t81 * t50 + (-cos(qJ(1)) * pkin(1) - t74 * t50 + t63 * t51) * qJD(1), 0, -t47 * t71 + t65 * t51 + t67, -t51 * t69 + (-t51 * t73 + t53 * t71) * pkin(4) + t67, t70, 0; t50 * qJD(5) - t81 * t51 + (-sin(qJ(1)) * pkin(1) + t74 * t51 + t63 * t50) * qJD(1), 0, t40 + t65 * t50 + (t47 - t64) * t70, t50 * t56 * t82 - t70 * t83 + t40, t71, 0; 0, 0, -t81, -t61, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:33
	% EndTime: 2019-10-10 01:23:33
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (532->65), mult. (430->88), div. (0->0), fcn. (321->12), ass. (0->57)
	t281 = qJ(3) + qJ(4);
	t275 = pkin(11) + t281;
	t271 = sin(t275);
	t284 = cos(qJ(6));
	t333 = r_i_i_C(1) * t284 + pkin(5);
	t297 = t333 * t271;
	t272 = cos(t275);
	t315 = qJD(6) * t271;
	t279 = qJD(3) + qJD(4);
	t282 = sin(qJ(6));
	t319 = t279 * t282;
	t336 = t272 * t319 + t284 * t315;
	t327 = pkin(9) + r_i_i_C(3);
	t335 = t327 * t272;
	t293 = -r_i_i_C(2) * t271 * t282 - t335;
	t276 = sin(t281);
	t283 = sin(qJ(3));
	t321 = pkin(3) * qJD(3);
	t312 = t283 * t321;
	t324 = pkin(4) * t279;
	t334 = -t276 * t324 + (-pkin(5) * t271 + t335) * t279 - t312;
	t304 = t282 * t315;
	t330 = r_i_i_C(1) * t304 + t336 * r_i_i_C(2);
	t298 = qJD(1) * t272 - qJD(6);
	t329 = t284 * t298;
	t314 = qJD(6) * t272;
	t299 = -qJD(1) + t314;
	t309 = t271 * t319;
	t328 = t299 * t284 - t309;
	t326 = pkin(4) * t276;
	t277 = cos(t281);
	t325 = pkin(4) * t277;
	t318 = t279 * t284;
	t280 = qJ(1) + pkin(10);
	t273 = sin(t280);
	t317 = qJD(1) * t273;
	t274 = cos(t280);
	t316 = qJD(1) * t274;
	t311 = t327 * t271;
	t296 = t330 * t274 + t317 * t297;
	t295 = t298 * t282;
	t294 = t330 * t273 - t293 * t316;
	t285 = cos(qJ(3));
	t292 = -t285 * pkin(3) - pkin(5) * t272 - pkin(2) - t311 - t325;
	t291 = -t297 - t326;
	t290 = -t272 * t333 - t311;
	t289 = t271 * t318 + t299 * t282;
	t288 = -t277 * t324 + t290 * t279 - t285 * t321;
	t287 = (t290 - t325) * t279;
	t286 = r_i_i_C(2) * t309 + (-r_i_i_C(1) * t282 - r_i_i_C(2) * t284) * t314 + (t291 + t335) * t279;
	t278 = -qJ(5) - pkin(8) - pkin(7);
	t270 = -t283 * pkin(3) - t326;
	t252 = t289 * t273 - t274 * t329;
	t251 = t328 * t273 + t274 * t295;
	t250 = t273 * t329 + t289 * t274;
	t249 = t273 * t295 - t328 * t274;
	t1 = [t252 * r_i_i_C(1) + t251 * r_i_i_C(2) + t274 * qJD(5) - t334 * t273 + (-cos(qJ(1)) * pkin(1) + t273 * t278 + t292 * t274) * qJD(1), 0, (-t270 + t293) * t317 + t288 * t274 + t296, (t293 + t326) * t317 + t274 * t287 + t296, t316, t249 * r_i_i_C(1) + t250 * r_i_i_C(2); -t250 * r_i_i_C(1) + t249 * r_i_i_C(2) + t273 * qJD(5) + t334 * t274 + (-sin(qJ(1)) * pkin(1) - t274 * t278 + t292 * t273) * qJD(1), 0, (t270 - t297) * t316 + t288 * t273 + t294, t273 * t287 + t291 * t316 + t294, t317, -t251 * r_i_i_C(1) + t252 * r_i_i_C(2); 0, 0, t286 - t312, t286, 0, (-t272 * t318 + t304) * r_i_i_C(2) - t336 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end