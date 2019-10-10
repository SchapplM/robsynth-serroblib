% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t106 = sin(qJ(2));
	t109 = cos(pkin(6)) * t106;
	t108 = qJD(2) * sin(pkin(7));
	t107 = cos(qJ(2));
	t104 = cos(pkin(12));
	t102 = sin(pkin(12));
	t1 = [0, 0, -(t102 * t109 - t104 * t107) * t108, 0, 0, 0; 0, 0, -(-t102 * t107 - t104 * t109) * t108, 0, 0, 0; 0, 0, sin(pkin(6)) * t106 * t108, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t141 = sin(qJ(2));
	t144 = cos(pkin(6)) * t141;
	t143 = qJD(2) * sin(pkin(7));
	t142 = cos(qJ(2));
	t139 = cos(pkin(12));
	t137 = sin(pkin(12));
	t1 = [0, 0, -(t137 * t144 - t139 * t142) * t143, 0, 0, 0; 0, 0, -(-t137 * t142 - t139 * t144) * t143, 0, 0, 0; 0, 0, sin(pkin(6)) * t141 * t143, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
	t183 = sin(pkin(7));
	t184 = sin(pkin(6));
	t202 = t184 * t183;
	t186 = cos(pkin(7));
	t188 = sin(qJ(3));
	t201 = t186 * t188;
	t187 = cos(pkin(6));
	t189 = sin(qJ(2));
	t200 = t187 * t189;
	t191 = cos(qJ(2));
	t199 = t187 * t191;
	t198 = t188 * t189;
	t190 = cos(qJ(3));
	t197 = t190 * t191;
	t196 = qJD(2) * t190;
	t182 = sin(pkin(12));
	t185 = cos(pkin(12));
	t195 = -t182 * t189 + t185 * t199;
	t194 = t182 * t191 + t185 * t200;
	t193 = -t182 * t199 - t185 * t189;
	t192 = t182 * t200 - t185 * t191;
	t181 = t192 * qJD(2);
	t180 = t194 * qJD(2);
	t1 = [0, 0, -t181 * t183, 0, t181 * t201 + t193 * t196 + (t192 * t188 + (t182 * t202 + t193 * t186) * t190) * qJD(3), 0; 0, 0, t180 * t183, 0, -t180 * t201 + t195 * t196 + (-t194 * t188 + (-t185 * t202 + t195 * t186) * t190) * qJD(3), 0; 0, 0, qJD(2) * t189 * t202, 0, t187 * t183 * qJD(3) * t190 + ((t186 * t197 - t198) * qJD(3) + (-t186 * t198 + t197) * qJD(2)) * t184, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:11
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (70->40), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->40)
	t263 = sin(qJ(3));
	t266 = cos(qJ(3));
	t256 = sin(pkin(12));
	t259 = cos(pkin(12));
	t267 = cos(qJ(2));
	t261 = cos(pkin(6));
	t264 = sin(qJ(2));
	t278 = t261 * t264;
	t269 = t256 * t278 - t259 * t267;
	t277 = t261 * t267;
	t254 = -t256 * t277 - t259 * t264;
	t260 = cos(pkin(7));
	t257 = sin(pkin(7));
	t258 = sin(pkin(6));
	t285 = t257 * t258;
	t270 = t254 * t260 + t256 * t285;
	t289 = t263 * t269 + t270 * t266;
	t253 = t256 * t267 + t259 * t278;
	t252 = -t256 * t264 + t259 * t277;
	t271 = -t252 * t260 + t259 * t285;
	t288 = t253 * t263 + t271 * t266;
	t262 = sin(qJ(5));
	t284 = t257 * t262;
	t283 = t257 * t264;
	t282 = t257 * t266;
	t281 = t258 * t260;
	t280 = t260 * t263;
	t279 = t260 * t266;
	t276 = t263 * t264;
	t275 = t263 * t267;
	t274 = t264 * t266;
	t273 = t266 * t267;
	t265 = cos(qJ(5));
	t272 = qJD(3) * t265;
	t268 = t260 * t273 - t276;
	t251 = t269 * qJD(2);
	t250 = t254 * qJD(2);
	t249 = t253 * qJD(2);
	t248 = t252 * qJD(2);
	t1 = [0, 0, -t251 * t257, 0, t289 * qJD(3) + t250 * t266 + t251 * t280, -t251 * t284 - (t250 * t263 - t251 * t279) * t265 + ((-t254 * t257 + t256 * t281) * t265 - t289 * t262) * qJD(5) - (t270 * t263 - t266 * t269) * t272; 0, 0, t249 * t257, 0, -t288 * qJD(3) + t248 * t266 - t249 * t280, t249 * t284 - (t248 * t263 + t249 * t279) * t265 + ((-t252 * t257 - t259 * t281) * t265 + t288 * t262) * qJD(5) - (t253 * t266 - t271 * t263) * t272; 0, 0, t258 * qJD(2) * t283, 0, t261 * qJD(3) * t282 + (t268 * qJD(3) + (-t260 * t276 + t273) * qJD(2)) * t258, (-t257 * t263 * t272 + (t260 * t265 - t262 * t282) * qJD(5)) * t261 + ((-t257 * t267 * t265 - t268 * t262) * qJD(5) - (t260 * t275 + t274) * t272 + (t262 * t283 - (t260 * t274 + t275) * t265) * qJD(2)) * t258;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end