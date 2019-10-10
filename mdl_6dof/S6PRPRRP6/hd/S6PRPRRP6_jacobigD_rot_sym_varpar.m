% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRP6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = cos(qJ(2));
	t83 = cos(pkin(6)) * t82;
	t81 = sin(qJ(2));
	t79 = cos(pkin(10));
	t78 = sin(pkin(10));
	t1 = [0, 0, 0, (-t78 * t83 - t79 * t81) * qJD(2), 0, 0; 0, 0, 0, (-t78 * t81 + t79 * t83) * qJD(2), 0, 0; 0, 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.05s
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
	t126 = sin(pkin(10));
	t128 = cos(pkin(10));
	t135 = -t126 * t131 + t128 * t137;
	t134 = t126 * t137 + t128 * t131;
	t130 = sin(qJ(4));
	t1 = [0, 0, 0, -t134 * qJD(2), (t126 * t140 + t134 * t130) * qJD(4) - (-t126 * t138 + t128 * t133) * t136, 0; 0, 0, 0, t135 * qJD(2), (-t128 * t140 - t135 * t130) * qJD(4) - (t126 * t133 + t128 * t138) * t136, 0; 0, 0, 0, qJD(2) * t139, -t127 * t131 * t136 + (t129 * t132 - t130 * t139) * qJD(4), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:57
	% EndTime: 2019-10-09 21:51:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (12->11), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t157 = sin(pkin(6));
	t162 = cos(qJ(4));
	t170 = t157 * t162;
	t163 = cos(qJ(2));
	t169 = t157 * t163;
	t159 = cos(pkin(6));
	t161 = sin(qJ(2));
	t168 = t159 * t161;
	t167 = t159 * t163;
	t166 = qJD(2) * t162;
	t156 = sin(pkin(10));
	t158 = cos(pkin(10));
	t165 = -t156 * t161 + t158 * t167;
	t164 = t156 * t167 + t158 * t161;
	t160 = sin(qJ(4));
	t1 = [0, 0, 0, -t164 * qJD(2), (t156 * t170 + t164 * t160) * qJD(4) - (-t156 * t168 + t158 * t163) * t166, 0; 0, 0, 0, t165 * qJD(2), (-t158 * t170 - t165 * t160) * qJD(4) - (t156 * t163 + t158 * t168) * t166, 0; 0, 0, 0, qJD(2) * t169, -t157 * t161 * t166 + (t159 * t162 - t160 * t169) * qJD(4), 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end