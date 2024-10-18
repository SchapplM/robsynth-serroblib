% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR13_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR13_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR13_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
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
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t48 = qJD(1) + qJD(2);
	t55 = pkin(2) * t48;
	t54 = pkin(1) * qJD(1);
	t49 = qJ(1) + qJ(2);
	t47 = qJ(3) + t49;
	t42 = sin(t47);
	t43 = cos(t47);
	t44 = qJD(3) + t48;
	t53 = (-r_i_i_C(1) * t43 + r_i_i_C(2) * t42) * t44;
	t52 = (-r_i_i_C(1) * t42 - r_i_i_C(2) * t43) * t44;
	t51 = -cos(t49) * t55 + t53;
	t50 = -sin(t49) * t55 + t52;
	t1 = [-cos(qJ(1)) * t54 + t51, t51, t53, 0, 0; -sin(qJ(1)) * t54 + t50, t50, t52, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (353->26), mult. (242->43), div. (0->0), fcn. (200->10), ass. (0->28)
	t179 = sin(pkin(5));
	t197 = t179 * (pkin(9) + r_i_i_C(3));
	t177 = qJD(1) + qJD(2);
	t195 = pkin(2) * t177;
	t194 = pkin(1) * qJD(1);
	t180 = cos(pkin(5));
	t181 = sin(qJ(4));
	t193 = t180 * t181;
	t182 = cos(qJ(4));
	t192 = t180 * t182;
	t178 = qJ(1) + qJ(2);
	t176 = qJ(3) + t178;
	t171 = sin(t176);
	t172 = cos(t176);
	t190 = t171 * t181 - t172 * t192;
	t189 = t171 * t182 + t172 * t193;
	t188 = t171 * t192 + t172 * t181;
	t187 = t171 * t193 - t172 * t182;
	t173 = qJD(3) + t177;
	t165 = t187 * qJD(4) + t190 * t173;
	t166 = t188 * qJD(4) + t189 * t173;
	t186 = -t166 * r_i_i_C(1) + t165 * r_i_i_C(2) + (-pkin(3) * t171 + t172 * t197) * t173;
	t185 = -sin(t178) * t195 + t186;
	t167 = t189 * qJD(4) + t188 * t173;
	t168 = t190 * qJD(4) + t187 * t173;
	t184 = t167 * r_i_i_C(2) + t168 * r_i_i_C(1) + (-pkin(3) * t172 - t171 * t197) * t173;
	t183 = -cos(t178) * t195 + t184;
	t1 = [-cos(qJ(1)) * t194 + t183, t183, t184, t165 * r_i_i_C(1) + t166 * r_i_i_C(2), 0; -sin(qJ(1)) * t194 + t185, t185, t186, -t167 * r_i_i_C(1) + t168 * r_i_i_C(2), 0; 0, 0, 0, (-r_i_i_C(1) * t181 - r_i_i_C(2) * t182) * t179 * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (831->52), mult. (451->75), div. (0->0), fcn. (314->16), ass. (0->52)
	t266 = qJD(4) + qJD(5);
	t299 = -t266 / 0.2e1;
	t268 = qJ(4) + qJ(5);
	t260 = pkin(5) - t268;
	t253 = sin(t260);
	t302 = t253 * t299;
	t259 = pkin(5) + t268;
	t297 = sin(t259);
	t244 = t297 / 0.2e1 - t253 / 0.2e1;
	t269 = qJ(1) + qJ(2);
	t265 = qJ(3) + t269;
	t255 = sin(t265);
	t256 = cos(t265);
	t267 = qJD(1) + qJD(2);
	t258 = qJD(3) + t267;
	t298 = cos(t260);
	t300 = cos(t259) / 0.2e1;
	t245 = t298 / 0.2e1 + t300;
	t263 = cos(t268);
	t283 = t245 * t266 + t258 * t263;
	t261 = sin(t268);
	t292 = t261 * t266;
	t301 = t283 * t255 + (t244 * t258 + t292) * t256;
	t296 = pkin(2) * t267;
	t295 = pkin(1) * qJD(1);
	t294 = pkin(4) * qJD(4);
	t293 = t255 * t258;
	t271 = cos(pkin(5));
	t272 = sin(qJ(4));
	t290 = t271 * t272;
	t273 = cos(qJ(4));
	t289 = t271 * t273;
	t280 = t245 * t258 + t263 * t266;
	t247 = t297 * t299;
	t282 = t258 * t261 - t247 + t302;
	t235 = t282 * t255 - t280 * t256;
	t288 = t235 * r_i_i_C(1) + t301 * r_i_i_C(2);
	t236 = t280 * t255 + t282 * t256;
	t278 = t244 * t293 + t255 * t292 - t283 * t256;
	t287 = -t236 * r_i_i_C(1) + t278 * r_i_i_C(2);
	t286 = (t266 * t300 + t298 * t299) * r_i_i_C(1) + (t247 + t302) * r_i_i_C(2);
	t270 = sin(pkin(5));
	t285 = r_i_i_C(3) * t258 * t270;
	t284 = t272 * t294;
	t279 = -t255 * t289 - t256 * t272;
	t246 = pkin(4) * t290 - t270 * (pkin(9) + pkin(10));
	t257 = t273 * pkin(4) + pkin(3);
	t277 = t236 * r_i_i_C(2) + t278 * r_i_i_C(1) + t246 * t293 + (-t257 * t258 - t289 * t294) * t256 + (-t285 + t284) * t255;
	t276 = t235 * r_i_i_C(2) - t301 * r_i_i_C(1) + t256 * t285 + (-t246 * t256 - t255 * t257) * t258 + t279 * t294;
	t275 = -cos(t269) * t296 + t277;
	t274 = -sin(t269) * t296 + t276;
	t1 = [-cos(qJ(1)) * t295 + t275, t275, t277, ((t255 * t272 - t256 * t289) * t258 + (t255 * t290 - t256 * t273) * qJD(4)) * pkin(4) + t288, t288; -sin(qJ(1)) * t295 + t274, t274, t276, (t279 * t258 + (-t255 * t273 - t256 * t290) * qJD(4)) * pkin(4) + t287, t287; 0, 0, 0, -t270 * t284 + t286, t286;];
	JaD_transl = t1;
end