% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = cos(qJ(2));
	t83 = cos(pkin(6)) * t82;
	t81 = sin(qJ(2));
	t79 = cos(pkin(11));
	t78 = sin(pkin(11));
	t1 = [0, 0, 0, (-t78 * t83 - t79 * t81) * qJD(2), 0, 0; 0, 0, 0, (-t78 * t81 + t79 * t83) * qJD(2), 0, 0; 0, 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (12->11), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t127 = sin(pkin(6));
	t132 = cos(qJ(4));
	t140 = t127 * t132;
	t133 = cos(qJ(2));
	t139 = t127 * t133;
	t129 = cos(pkin(6));
	t131 = sin(qJ(2));
	t138 = t129 * t131;
	t137 = t129 * t133;
	t136 = qJD(2) * t132;
	t126 = sin(pkin(11));
	t128 = cos(pkin(11));
	t135 = -t126 * t131 + t128 * t137;
	t134 = t126 * t137 + t128 * t131;
	t130 = sin(qJ(4));
	t1 = [0, 0, 0, -t134 * qJD(2), (t126 * t140 + t134 * t130) * qJD(4) - (-t126 * t138 + t128 * t133) * t136, 0; 0, 0, 0, t135 * qJD(2), (-t128 * t140 - t135 * t130) * qJD(4) - (t126 * t133 + t128 * t138) * t136, 0; 0, 0, 0, qJD(2) * t139, -t127 * t131 * t136 + (t129 * t132 - t130 * t139) * qJD(4), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (22->11), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
	t155 = sin(pkin(6));
	t160 = cos(qJ(4));
	t168 = t155 * t160;
	t161 = cos(qJ(2));
	t167 = t155 * t161;
	t157 = cos(pkin(6));
	t159 = sin(qJ(2));
	t166 = t157 * t159;
	t165 = t157 * t161;
	t164 = qJD(2) * t160;
	t154 = sin(pkin(11));
	t156 = cos(pkin(11));
	t163 = -t154 * t159 + t156 * t165;
	t162 = t154 * t165 + t156 * t159;
	t158 = sin(qJ(4));
	t153 = -t155 * t159 * t164 + (t157 * t160 - t158 * t167) * qJD(4);
	t152 = (-t156 * t168 - t163 * t158) * qJD(4) - (t154 * t161 + t156 * t166) * t164;
	t151 = (t154 * t168 + t162 * t158) * qJD(4) - (-t154 * t166 + t156 * t161) * t164;
	t1 = [0, 0, 0, -t162 * qJD(2), t151, t151; 0, 0, 0, t163 * qJD(2), t152, t152; 0, 0, 0, qJD(2) * t167, t153, t153;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end