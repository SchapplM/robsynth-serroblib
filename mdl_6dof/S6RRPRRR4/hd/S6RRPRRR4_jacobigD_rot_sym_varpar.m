% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR4_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t79 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t79, 0, 0, 0, 0; 0, sin(qJ(1)) * t79, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t133 = sin(pkin(12));
	t135 = cos(pkin(12));
	t137 = sin(qJ(2));
	t139 = cos(qJ(2));
	t142 = t139 * t133 + t137 * t135;
	t144 = qJD(2) * t142;
	t134 = sin(pkin(6));
	t143 = qJD(1) * t134;
	t141 = t133 * t137 - t135 * t139;
	t140 = cos(qJ(1));
	t138 = sin(qJ(1));
	t136 = cos(pkin(6));
	t131 = t141 * qJD(2);
	t130 = t141 * t136;
	t129 = t136 * t144;
	t1 = [0, t140 * t143, 0, -t138 * t129 - t140 * t131 + (-t130 * t140 - t138 * t142) * qJD(1), 0, 0; 0, t138 * t143, 0, t140 * t129 - t138 * t131 + (-t130 * t138 + t140 * t142) * qJD(1), 0, 0; 0, 0, 0, t134 * t144, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (30->8), mult. (104->22), div. (0->0), fcn. (112->8), ass. (0->19)
	t158 = sin(pkin(12));
	t160 = cos(pkin(12));
	t162 = sin(qJ(2));
	t164 = cos(qJ(2));
	t167 = t164 * t158 + t162 * t160;
	t169 = qJD(2) * t167;
	t159 = sin(pkin(6));
	t168 = qJD(1) * t159;
	t166 = t158 * t162 - t160 * t164;
	t165 = cos(qJ(1));
	t163 = sin(qJ(1));
	t161 = cos(pkin(6));
	t156 = t166 * qJD(2);
	t155 = t166 * t161;
	t154 = t161 * t169;
	t153 = t159 * t169;
	t152 = t165 * t154 - t163 * t156 + (-t155 * t163 + t165 * t167) * qJD(1);
	t151 = -t163 * t154 - t165 * t156 + (-t155 * t165 - t163 * t167) * qJD(1);
	t1 = [0, t165 * t168, 0, t151, t151, 0; 0, t163 * t168, 0, t152, t152, 0; 0, 0, 0, t153, t153, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (77->23), mult. (197->50), div. (0->0), fcn. (216->10), ass. (0->29)
	t242 = sin(pkin(12));
	t244 = cos(pkin(12));
	t246 = sin(qJ(2));
	t248 = cos(qJ(2));
	t251 = t246 * t242 - t248 * t244;
	t259 = t251 * qJD(2);
	t252 = t248 * t242 + t246 * t244;
	t235 = t252 * qJD(2);
	t241 = qJ(4) + qJ(5);
	t239 = cos(t241);
	t240 = qJD(4) + qJD(5);
	t258 = t239 * t240;
	t243 = sin(pkin(6));
	t257 = t243 * t239;
	t256 = qJD(1) * t243;
	t245 = cos(pkin(6));
	t255 = t240 * t243 + t245 * t259;
	t233 = t252 * t245;
	t247 = sin(qJ(1));
	t249 = cos(qJ(1));
	t254 = t249 * t233 - t247 * t251;
	t253 = -t247 * t233 - t249 * t251;
	t238 = sin(t241);
	t232 = t251 * t245;
	t230 = t245 * t235;
	t229 = t243 * t235;
	t228 = t249 * t230 - t247 * t259 + (-t232 * t247 + t249 * t252) * qJD(1);
	t227 = -t247 * t230 - t249 * t259 + (-t232 * t249 - t247 * t252) * qJD(1);
	t1 = [0, t249 * t256, 0, t227, t227, t253 * t258 + (-t249 * t235 + t255 * t247) * t238 + (-t254 * t238 - t249 * t257) * qJD(1); 0, t247 * t256, 0, t228, t228, t254 * t258 + (-t247 * t235 - t255 * t249) * t238 + (t253 * t238 - t247 * t257) * qJD(1); 0, 0, 0, t229, t229, t245 * t240 * t238 + (-t238 * t259 + t252 * t258) * t243;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end