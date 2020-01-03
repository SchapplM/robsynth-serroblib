% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:15:13
	% EndTime: 2019-12-31 21:15:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:15:14
	% EndTime: 2019-12-31 21:15:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:15:13
	% EndTime: 2019-12-31 21:15:13
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:15:14
	% EndTime: 2019-12-31 21:15:14
	% DurationCPUTime: 0.10s
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
	t58 = r_i_i_C(3) + pkin(7) + pkin(6);
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
	t1 = [t63 * t41 + (-t41 * t58 + t43 * t48) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0; -t63 * t43 + (t48 * t41 + t43 * t58) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0; 0, -t63, -t47, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:15:14
	% EndTime: 2019-12-31 21:15:14
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (149->33), mult. (146->42), div. (0->0), fcn. (95->8), ass. (0->32)
	t52 = sin(qJ(2));
	t50 = qJD(2) + qJD(3);
	t51 = qJ(2) + qJ(3);
	t46 = pkin(9) + t51;
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
	t69 = r_i_i_C(3) + qJ(4) + pkin(7) + pkin(6);
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
	t1 = [t55 * qJD(4) + t76 * t53 + (-t69 * t53 + t58 * t55) * qJD(1), -t43 * t66 + t60 * t55 + t62, -t55 * t64 + (t47 * t66 - t55 * t68) * pkin(3) + t62, t65, 0; t53 * qJD(4) - t76 * t55 + (t58 * t53 + t69 * t55) * qJD(1), t38 + t60 * t53 + (t43 - t59) * t65, t50 * t53 * t77 - t65 * t78 + t38, t66, 0; 0, -t76, -t56, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:15:15
	% EndTime: 2019-12-31 21:15:15
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (412->63), mult. (426->88), div. (0->0), fcn. (319->10), ass. (0->56)
	t271 = qJ(2) + qJ(3);
	t266 = pkin(9) + t271;
	t264 = sin(t266);
	t275 = cos(qJ(5));
	t322 = r_i_i_C(1) * t275 + pkin(4);
	t289 = t322 * t264;
	t272 = sin(qJ(5));
	t306 = qJD(5) * t264;
	t265 = cos(t266);
	t270 = qJD(2) + qJD(3);
	t311 = t265 * t270;
	t325 = t272 * t311 + t275 * t306;
	t319 = pkin(8) + r_i_i_C(3);
	t301 = t319 * t265;
	t304 = r_i_i_C(2) * t264 * t272;
	t324 = t301 + t304;
	t267 = sin(t271);
	t273 = sin(qJ(2));
	t313 = pkin(2) * qJD(2);
	t303 = t273 * t313;
	t316 = pkin(3) * t270;
	t323 = -t267 * t316 + (-pkin(4) * t264 + t301) * t270 - t303;
	t296 = t272 * t306;
	t320 = r_i_i_C(1) * t296 + t325 * r_i_i_C(2);
	t318 = pkin(3) * t267;
	t268 = cos(t271);
	t317 = pkin(3) * t268;
	t312 = t264 * t270;
	t310 = t270 * t275;
	t277 = cos(qJ(1));
	t309 = t275 * t277;
	t274 = sin(qJ(1));
	t308 = qJD(1) * t274;
	t307 = qJD(1) * t277;
	t305 = qJD(5) * t265;
	t302 = t319 * t264;
	t291 = -qJD(1) + t305;
	t290 = qJD(1) * t265 - qJD(5);
	t288 = t320 * t277 + t308 * t289;
	t287 = t291 * t272;
	t286 = t320 * t274 + t324 * t307;
	t276 = cos(qJ(2));
	t284 = -t276 * pkin(2) - pkin(4) * t265 - pkin(1) - t302 - t317;
	t283 = -t289 - t318;
	t282 = -t265 * t322 - t302;
	t281 = t290 * t274 + t277 * t312;
	t280 = -t268 * t316 + t282 * t270 - t276 * t313;
	t279 = t270 * (t282 - t317);
	t278 = (-r_i_i_C(1) * t272 - r_i_i_C(2) * t275) * t305 + t319 * t311 + (t304 + t283) * t270;
	t269 = -qJ(4) - pkin(7) - pkin(6);
	t263 = -t273 * pkin(2) - t318;
	t245 = -t290 * t309 + (t264 * t310 + t287) * t274;
	t244 = t291 * t275 * t274 + (-t274 * t312 + t290 * t277) * t272;
	t243 = t281 * t275 + t277 * t287;
	t242 = t281 * t272 - t291 * t309;
	t1 = [t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t277 * qJD(4) - t323 * t274 + (t269 * t274 + t284 * t277) * qJD(1), (-t263 - t324) * t308 + t280 * t277 + t288, (-t324 + t318) * t308 + t277 * t279 + t288, t307, t242 * r_i_i_C(1) + t243 * r_i_i_C(2); -t243 * r_i_i_C(1) + t242 * r_i_i_C(2) + t274 * qJD(4) + t323 * t277 + (-t269 * t277 + t284 * t274) * qJD(1), (t263 - t289) * t307 + t280 * t274 + t286, t274 * t279 + t283 * t307 + t286, t308, -t244 * r_i_i_C(1) + t245 * r_i_i_C(2); 0, t278 - t303, t278, 0, (-t265 * t310 + t296) * r_i_i_C(2) - t325 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end