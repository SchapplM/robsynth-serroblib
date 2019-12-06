% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiaD_transl_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:12
	% EndTime: 2019-12-05 18:54:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:12
	% EndTime: 2019-12-05 18:54:12
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
	% StartTime: 2019-12-05 18:54:13
	% EndTime: 2019-12-05 18:54:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:13
	% EndTime: 2019-12-05 18:54:13
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (69->16), mult. (88->32), div. (0->0), fcn. (56->6), ass. (0->18)
	t36 = qJ(1) + qJ(2);
	t33 = sin(t36);
	t34 = cos(t36);
	t38 = cos(qJ(3));
	t47 = qJD(3) * t38;
	t35 = qJD(1) + qJD(2);
	t37 = sin(qJ(3));
	t51 = t35 * t37;
	t53 = t33 * t51 - t34 * t47;
	t52 = t33 * t47 + t34 * t51;
	t50 = t35 * t38;
	t49 = pkin(1) * qJD(1);
	t48 = qJD(3) * t37;
	t44 = t33 * t48;
	t41 = t33 * t50 + t34 * t48;
	t40 = r_i_i_C(1) * t44 + (-r_i_i_C(1) * t34 * t38 - r_i_i_C(3) * t33) * t35 + t52 * r_i_i_C(2);
	t39 = t35 * t34 * r_i_i_C(3) - t41 * r_i_i_C(1) + t53 * r_i_i_C(2);
	t1 = [-cos(qJ(1)) * t49 + t40, t40, t53 * r_i_i_C(1) + t41 * r_i_i_C(2), 0, 0; -sin(qJ(1)) * t49 + t39, t39, (-t34 * t50 + t44) * r_i_i_C(2) - t52 * r_i_i_C(1), 0, 0; 0, 0, (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:13
	% EndTime: 2019-12-05 18:54:13
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (175->28), mult. (154->44), div. (0->0), fcn. (99->8), ass. (0->33)
	t57 = qJ(3) + qJ(4);
	t51 = sin(t57);
	t58 = qJ(1) + qJ(2);
	t52 = sin(t58);
	t54 = cos(t58);
	t56 = qJD(1) + qJD(2);
	t77 = t54 * t56;
	t53 = cos(t57);
	t55 = qJD(3) + qJD(4);
	t78 = t53 * t55;
	t83 = t51 * t77 + t52 * t78;
	t59 = sin(qJ(3));
	t82 = pkin(2) * t59;
	t81 = r_i_i_C(2) * t53;
	t80 = t51 * t55;
	t79 = t52 * t56;
	t76 = pkin(1) * qJD(1);
	t60 = cos(qJ(3));
	t75 = qJD(3) * t60;
	t74 = r_i_i_C(1) * t78;
	t73 = t56 * t81;
	t72 = t52 * t80;
	t71 = t51 * t79;
	t68 = qJD(3) * t82;
	t67 = -pkin(2) * t60 - r_i_i_C(1) * t53;
	t66 = -r_i_i_C(1) * t51 - t81;
	t65 = t66 * t55;
	t64 = r_i_i_C(1) * t71 + t52 * t73 + (r_i_i_C(2) * t80 - t74) * t54;
	t63 = t65 - t68;
	t62 = r_i_i_C(1) * t72 + t52 * t68 + (-r_i_i_C(3) * t52 + t67 * t54) * t56 + t83 * r_i_i_C(2);
	t61 = r_i_i_C(2) * t71 + r_i_i_C(3) * t77 + t63 * t54 + t67 * t79;
	t41 = r_i_i_C(2) * t72;
	t1 = [-cos(qJ(1)) * t76 + t62, t62, (-t54 * t75 + t59 * t79) * pkin(2) + t64, t64, 0; -sin(qJ(1)) * t76 + t61, t61, t41 + (-pkin(2) * t75 - t74) * t52 + (t66 - t82) * t77, -t83 * r_i_i_C(1) - t54 * t73 + t41, 0; 0, 0, t63, t65, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:14
	% EndTime: 2019-12-05 18:54:14
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (400->53), mult. (394->92), div. (0->0), fcn. (313->10), ass. (0->56)
	t274 = qJ(3) + qJ(4);
	t268 = sin(t274);
	t270 = cos(t274);
	t278 = cos(qJ(5));
	t304 = qJD(5) * t278;
	t272 = qJD(3) + qJD(4);
	t276 = sin(qJ(5));
	t309 = t272 * t276;
	t322 = t268 * t304 + t270 * t309;
	t305 = qJD(5) * t276;
	t296 = t268 * t305;
	t321 = r_i_i_C(1) * t296 + t322 * r_i_i_C(2);
	t273 = qJD(1) + qJD(2);
	t311 = t270 * t273;
	t292 = -qJD(5) + t311;
	t320 = t278 * t292;
	t291 = qJD(5) * t270 - t273;
	t300 = t268 * t309;
	t319 = t291 * t278 - t300;
	t277 = sin(qJ(3));
	t318 = pkin(2) * t277;
	t317 = r_i_i_C(1) * t278;
	t316 = r_i_i_C(3) * t270;
	t315 = pkin(1) * qJD(1);
	t314 = t268 * t272;
	t275 = qJ(1) + qJ(2);
	t269 = sin(t275);
	t313 = t269 * t273;
	t312 = t270 * t272;
	t271 = cos(t275);
	t310 = t271 * t273;
	t308 = t272 * t278;
	t279 = cos(qJ(3));
	t307 = t273 * t279;
	t306 = qJD(3) * t277;
	t303 = r_i_i_C(2) * t268 * t276;
	t267 = r_i_i_C(3) * t312;
	t302 = pkin(2) * t306;
	t301 = t268 * t313;
	t299 = t268 * t308;
	t297 = t270 * t308;
	t294 = t273 * t303;
	t288 = t321 * t271 + t301 * t317;
	t287 = t321 * t269 + t271 * t294 + t310 * t316;
	t286 = t292 * t276;
	t285 = -t268 * t310 - t269 * t312;
	t284 = t291 * t276 + t299;
	t283 = -pkin(2) * qJD(3) * t279 + (-r_i_i_C(3) * t268 - t270 * t317) * t272;
	t282 = t267 + (-t270 * t305 - t299) * r_i_i_C(1) + (-t270 * t304 + t300) * r_i_i_C(2);
	t253 = t319 * t269 + t271 * t286;
	t254 = t284 * t269 - t271 * t320;
	t281 = -t271 * pkin(2) * t307 + t254 * r_i_i_C(1) + t253 * r_i_i_C(2) + t285 * r_i_i_C(3) + t269 * t302;
	t251 = t269 * t286 - t319 * t271;
	t252 = t269 * t320 + t284 * t271;
	t280 = -r_i_i_C(3) * t301 + t251 * r_i_i_C(2) - t252 * r_i_i_C(1) + t271 * t267 + (-t269 * t307 - t271 * t306) * pkin(2);
	t1 = [-cos(qJ(1)) * t315 + t281, t281, t283 * t271 + (-t303 - t316 + t318) * t313 + t288, -t271 * r_i_i_C(1) * t297 - t269 * t294 + (-t269 * t311 - t271 * t314) * r_i_i_C(3) + t288, t251 * r_i_i_C(1) + t252 * r_i_i_C(2); -sin(qJ(1)) * t315 + t280, t280, (-t268 * t317 - t318) * t310 + t283 * t269 + t287, -t269 * r_i_i_C(3) * t314 + t285 * t317 + t287, -t253 * r_i_i_C(1) + t254 * r_i_i_C(2); 0, 0, t282 - t302, t282, (t296 - t297) * r_i_i_C(2) - t322 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end