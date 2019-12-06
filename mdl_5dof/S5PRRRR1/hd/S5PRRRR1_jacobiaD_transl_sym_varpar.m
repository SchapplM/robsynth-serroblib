% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR1
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiaD_transl_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:49
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:49
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:49
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(2));
	t5 = sin(qJ(2));
	t1 = [0, (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0; 0, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:49
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->12), mult. (56->29), div. (0->0), fcn. (36->4), ass. (0->12)
	t66 = sin(qJ(2));
	t75 = qJD(2) * t66;
	t68 = cos(qJ(2));
	t74 = qJD(2) * t68;
	t73 = qJD(3) * t66;
	t72 = qJD(3) * t68;
	t65 = sin(qJ(3));
	t67 = cos(qJ(3));
	t71 = -r_i_i_C(1) * t67 + r_i_i_C(2) * t65;
	t70 = r_i_i_C(1) * t65 + r_i_i_C(2) * t67;
	t69 = t70 * qJD(3);
	t1 = [0, t66 * t69 + (-r_i_i_C(3) * t66 + t71 * t68) * qJD(2), (t65 * t72 + t67 * t75) * r_i_i_C(2) + (t65 * t75 - t67 * t72) * r_i_i_C(1), 0, 0; 0, 0, t69, 0, 0; 0, -t70 * t72 + (r_i_i_C(3) * t68 + t71 * t66) * qJD(2), (t65 * t73 - t67 * t74) * r_i_i_C(2) + (-t65 * t74 - t67 * t73) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:49
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (75->24), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
	t93 = sin(qJ(3));
	t112 = pkin(2) * t93;
	t100 = qJD(3) * t112;
	t91 = qJD(3) + qJD(4);
	t92 = qJ(3) + qJ(4);
	t90 = cos(t92);
	t110 = r_i_i_C(2) * t90;
	t89 = sin(t92);
	t111 = r_i_i_C(1) * t89;
	t99 = t110 + t111;
	t113 = t99 * t91 + t100;
	t109 = t89 * t91;
	t108 = t90 * t91;
	t107 = r_i_i_C(1) * t109 + r_i_i_C(2) * t108;
	t94 = sin(qJ(2));
	t106 = qJD(2) * t94;
	t96 = cos(qJ(2));
	t105 = qJD(2) * t96;
	t95 = cos(qJ(3));
	t104 = qJD(3) * t95;
	t103 = r_i_i_C(1) * t108;
	t102 = r_i_i_C(2) * t109;
	t101 = qJD(2) * t110;
	t98 = -pkin(2) * t95 - r_i_i_C(1) * t90 + r_i_i_C(2) * t89;
	t97 = t94 * t101 + t106 * t111 + (t102 - t103) * t96;
	t83 = t94 * t102;
	t1 = [0, t113 * t94 + (-r_i_i_C(3) * t94 + t98 * t96) * qJD(2), (-t96 * t104 + t93 * t106) * pkin(2) + t97, t97, 0; 0, 0, t100 + t107, t107, 0; 0, -t113 * t96 + (r_i_i_C(3) * t96 + t98 * t94) * qJD(2), t83 + (-pkin(2) * t104 - t103) * t94 + (-t99 - t112) * t105, -t96 * t101 + t83 + (-t89 * t105 - t94 * t108) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:50
	% EndTime: 2019-12-05 17:03:50
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (192->49), mult. (310->87), div. (0->0), fcn. (243->8), ass. (0->47)
	t252 = qJ(3) + qJ(4);
	t249 = sin(t252);
	t250 = cos(t252);
	t251 = qJD(3) + qJD(4);
	t253 = sin(qJ(5));
	t256 = cos(qJ(5));
	t282 = qJD(5) * t256;
	t295 = t250 * t251 * t253 + t249 * t282;
	t283 = qJD(5) * t253;
	t274 = t249 * t283;
	t294 = r_i_i_C(1) * t274 + t295 * r_i_i_C(2);
	t254 = sin(qJ(3));
	t293 = pkin(2) * t254;
	t292 = r_i_i_C(1) * t256;
	t291 = r_i_i_C(3) * t250;
	t255 = sin(qJ(2));
	t290 = t251 * t255;
	t289 = t251 * t256;
	t258 = cos(qJ(2));
	t288 = t251 * t258;
	t287 = t256 * t258;
	t286 = qJD(2) * t255;
	t285 = qJD(2) * t258;
	t284 = qJD(3) * t254;
	t281 = r_i_i_C(2) * t249 * t253;
	t280 = t249 * t290;
	t279 = t249 * t289;
	t278 = t249 * t288;
	t276 = t250 * t288;
	t275 = t249 * t286;
	t270 = qJD(2) * t281;
	t268 = qJD(5) * t250 - qJD(2);
	t267 = qJD(2) * t250 - qJD(5);
	t266 = t294 * t258 + t275 * t292;
	t265 = t294 * t255 + t258 * t270 + t285 * t291;
	t264 = t268 * t253;
	t263 = -t281 - t291;
	t262 = -t249 * t285 - t250 * t290;
	t261 = t255 * t267 + t278;
	t260 = t250 * r_i_i_C(2) * t282 + t251 * t263 + (t250 * t283 + t279) * r_i_i_C(1);
	t257 = cos(qJ(3));
	t259 = -pkin(2) * qJD(3) * t257 + (-r_i_i_C(3) * t249 - t250 * t292) * t251;
	t236 = -t267 * t287 + (t264 + t279) * t255;
	t235 = t268 * t256 * t255 + (t267 * t258 - t280) * t253;
	t234 = t256 * t261 + t258 * t264;
	t233 = t253 * t261 - t268 * t287;
	t1 = [0, t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t262 * r_i_i_C(3) + (t255 * t284 - t257 * t285) * pkin(2), t259 * t258 + (t263 + t293) * t286 + t266, -t276 * t292 - t255 * t270 + (-t250 * t286 - t278) * r_i_i_C(3) + t266, t233 * r_i_i_C(1) + t234 * r_i_i_C(2); 0, 0, pkin(2) * t284 + t260, t260, (t250 * t289 - t274) * r_i_i_C(2) + t295 * r_i_i_C(1); 0, -t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + (-t275 + t276) * r_i_i_C(3) + (-t257 * t286 - t258 * t284) * pkin(2), (-t249 * t292 - t293) * t285 + t259 * t255 + t265, -r_i_i_C(3) * t280 + t262 * t292 + t265, -t235 * r_i_i_C(1) + t236 * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end