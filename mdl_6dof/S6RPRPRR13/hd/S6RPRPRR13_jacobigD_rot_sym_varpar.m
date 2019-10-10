% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR13_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t127 = sin(pkin(6)) * cos(pkin(7));
	t126 = cos(pkin(6)) * cos(pkin(12));
	t125 = cos(qJ(1));
	t124 = sin(qJ(1));
	t119 = sin(pkin(7));
	t118 = sin(pkin(12));
	t1 = [0, 0, (-(t118 * t124 - t125 * t126) * t119 + t125 * t127) * qJD(1), 0, 0, 0; 0, 0, (-(-t118 * t125 - t124 * t126) * t119 + t124 * t127) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t151 = sin(pkin(6)) * cos(pkin(7));
	t150 = cos(pkin(6)) * cos(pkin(12));
	t149 = cos(qJ(1));
	t148 = sin(qJ(1));
	t143 = sin(pkin(7));
	t142 = sin(pkin(12));
	t1 = [0, 0, (-(t142 * t148 - t149 * t150) * t143 + t149 * t151) * qJD(1), 0, 0, 0; 0, 0, (-(-t142 * t149 - t148 * t150) * t143 + t148 * t151) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (24->18), mult. (92->43), div. (0->0), fcn. (96->10), ass. (0->25)
	t185 = sin(pkin(7));
	t186 = sin(pkin(6));
	t205 = t186 * t185;
	t188 = cos(pkin(7));
	t190 = sin(qJ(3));
	t204 = t188 * t190;
	t189 = cos(pkin(6));
	t193 = cos(qJ(1));
	t203 = t189 * t193;
	t184 = sin(pkin(12));
	t191 = sin(qJ(1));
	t202 = t191 * t184;
	t187 = cos(pkin(12));
	t201 = t191 * t187;
	t200 = t191 * t205;
	t199 = t193 * t205;
	t198 = qJD(1) * t186 * t188;
	t197 = t187 * t203 - t202;
	t196 = -t184 * t193 - t189 * t201;
	t195 = t184 * t203 + t201;
	t194 = t187 * t193 - t189 * t202;
	t192 = cos(qJ(3));
	t183 = t196 * qJD(1);
	t182 = t197 * qJD(1);
	t1 = [0, 0, t182 * t185 + t193 * t198, 0, -t182 * t204 + (-t194 * t190 + (t196 * t188 + t200) * t192) * qJD(3) + (t190 * t199 - t195 * t192) * qJD(1), 0; 0, 0, -t183 * t185 + t191 * t198, 0, t183 * t204 + (-t195 * t190 + (t197 * t188 - t199) * t192) * qJD(3) + (t190 * t200 + t194 * t192) * qJD(1), 0; 0, 0, 0, 0, (t185 * t189 * t192 + (t187 * t188 * t192 - t184 * t190) * t186) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:46
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
	t258 = sin(pkin(12));
	t260 = sin(pkin(6));
	t261 = cos(pkin(12));
	t265 = sin(qJ(3));
	t268 = cos(qJ(3));
	t262 = cos(pkin(7));
	t284 = t262 * t268;
	t259 = sin(pkin(7));
	t263 = cos(pkin(6));
	t287 = t259 * t263;
	t292 = (-t258 * t265 + t261 * t284) * t260 + t268 * t287;
	t269 = cos(qJ(1));
	t280 = t269 * t261;
	t266 = sin(qJ(1));
	t283 = t266 * t258;
	t257 = -t263 * t283 + t280;
	t281 = t269 * t258;
	t282 = t266 * t261;
	t256 = -t263 * t282 - t281;
	t286 = t260 * t266;
	t271 = t256 * t262 + t259 * t286;
	t291 = -t257 * t265 + t271 * t268;
	t255 = t263 * t281 + t282;
	t254 = t263 * t280 - t283;
	t285 = t260 * t269;
	t272 = -t254 * t262 + t259 * t285;
	t290 = t255 * t265 + t272 * t268;
	t279 = qJD(1) * t260;
	t267 = cos(qJ(5));
	t278 = qJD(3) * t267;
	t276 = t266 * t279;
	t275 = t269 * t279;
	t274 = t259 * t276;
	t273 = t259 * t275;
	t264 = sin(qJ(5));
	t253 = t257 * qJD(1);
	t252 = t256 * qJD(1);
	t251 = t255 * qJD(1);
	t250 = t254 * qJD(1);
	t249 = -t252 * t259 + t262 * t276;
	t248 = t250 * t259 + t262 * t275;
	t1 = [0, 0, t248, 0, -t251 * t268 + (-t250 * t262 + t273) * t265 + t291 * qJD(3), t248 * t264 - (t250 * t284 - t251 * t265 - t268 * t273) * t267 + ((-t256 * t259 + t262 * t286) * t267 - t291 * t264) * qJD(5) - (t257 * t268 + t271 * t265) * t278; 0, 0, t249, 0, t253 * t268 + (t252 * t262 + t274) * t265 - t290 * qJD(3), t249 * t264 - (-t252 * t284 + t253 * t265 - t268 * t274) * t267 + ((-t254 * t259 - t262 * t285) * t267 + t290 * t264) * qJD(5) - (t255 * t268 - t272 * t265) * t278; 0, 0, 0, 0, t292 * qJD(3), ((-t260 * t261 * t259 + t263 * t262) * t267 - t292 * t264) * qJD(5) - (t265 * t287 + (t261 * t262 * t265 + t258 * t268) * t260) * t278;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end