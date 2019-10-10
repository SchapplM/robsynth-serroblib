% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->10), mult. (40->28), div. (0->0), fcn. (44->10), ass. (0->14)
	t136 = sin(pkin(12));
	t142 = cos(pkin(6));
	t147 = t136 * t142;
	t137 = sin(pkin(7));
	t138 = sin(pkin(6));
	t146 = t138 * t137;
	t140 = cos(pkin(12));
	t145 = t140 * t142;
	t144 = cos(qJ(3));
	t143 = sin(qJ(3));
	t141 = cos(pkin(7));
	t139 = cos(pkin(13));
	t135 = sin(pkin(13));
	t1 = [0, 0, 0, ((-t135 * t147 + t140 * t139) * t144 + ((-t140 * t135 - t139 * t147) * t141 + t136 * t146) * t143) * qJD(3), 0, 0; 0, 0, 0, ((t135 * t145 + t136 * t139) * t144 + ((-t136 * t135 + t139 * t145) * t141 - t140 * t146) * t143) * qJD(3), 0, 0; 0, 0, 0, (t137 * t142 * t143 + (t139 * t141 * t143 + t135 * t144) * t138) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (20->10), mult. (80->28), div. (0->0), fcn. (88->10), ass. (0->17)
	t165 = sin(pkin(12));
	t171 = cos(pkin(6));
	t176 = t165 * t171;
	t166 = sin(pkin(7));
	t167 = sin(pkin(6));
	t175 = t167 * t166;
	t169 = cos(pkin(12));
	t174 = t169 * t171;
	t173 = cos(qJ(3));
	t172 = sin(qJ(3));
	t170 = cos(pkin(7));
	t168 = cos(pkin(13));
	t164 = sin(pkin(13));
	t163 = (t166 * t171 * t172 + (t168 * t170 * t172 + t164 * t173) * t167) * qJD(3);
	t162 = ((-t164 * t176 + t169 * t168) * t173 + ((-t169 * t164 - t168 * t176) * t170 + t165 * t175) * t172) * qJD(3);
	t161 = ((t164 * t174 + t165 * t168) * t173 + ((-t165 * t164 + t168 * t174) * t170 - t169 * t175) * t172) * qJD(3);
	t1 = [0, 0, 0, t162, t162, 0; 0, 0, 0, t161, t161, 0; 0, 0, 0, t163, t163, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (66->25), mult. (181->56), div. (0->0), fcn. (208->12), ass. (0->34)
	t252 = sin(pkin(12));
	t258 = cos(pkin(6));
	t272 = t252 * t258;
	t253 = sin(pkin(7));
	t254 = sin(pkin(6));
	t271 = t253 * t254;
	t270 = t253 * t258;
	t257 = cos(pkin(7));
	t269 = t254 * t257;
	t255 = cos(pkin(13));
	t268 = t255 * t257;
	t256 = cos(pkin(12));
	t267 = t256 * t258;
	t250 = qJ(4) + qJ(5);
	t247 = sin(t250);
	t266 = qJD(3) * t247;
	t251 = sin(pkin(13));
	t243 = -t252 * t251 + t255 * t267;
	t265 = t243 * t257 - t256 * t271;
	t245 = -t256 * t251 - t255 * t272;
	t264 = t245 * t257 + t252 * t271;
	t244 = t251 * t267 + t252 * t255;
	t259 = sin(qJ(3));
	t260 = cos(qJ(3));
	t263 = t244 * t260 + t265 * t259;
	t246 = -t251 * t272 + t256 * t255;
	t262 = t246 * t260 + t264 * t259;
	t261 = t259 * t270 + (t251 * t260 + t259 * t268) * t254;
	t249 = qJD(4) + qJD(5);
	t248 = cos(t250);
	t242 = t261 * qJD(3);
	t241 = t262 * qJD(3);
	t240 = t263 * qJD(3);
	t1 = [0, 0, 0, t241, t241, (t262 * t248 + (-t245 * t253 + t252 * t269) * t247) * t249 + (-t246 * t259 + t264 * t260) * t266; 0, 0, 0, t240, t240, (t263 * t248 + (-t243 * t253 - t256 * t269) * t247) * t249 + (-t244 * t259 + t265 * t260) * t266; 0, 0, 0, t242, t242, (t261 * t248 + (-t255 * t271 + t258 * t257) * t247) * t249 + (t260 * t270 + (-t251 * t259 + t260 * t268) * t254) * t266;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end