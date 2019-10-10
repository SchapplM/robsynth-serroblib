% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRP4_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP4_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(10));
	t79 = sin(pkin(10));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:56
	% EndTime: 2019-10-09 22:21:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t97 = sin(qJ(2));
	t99 = cos(pkin(6)) * t97;
	t98 = cos(qJ(2));
	t95 = cos(pkin(10));
	t94 = sin(pkin(10));
	t1 = [0, 0, (-t94 * t99 + t95 * t98) * qJD(2), 0, 0, 0; 0, 0, (t94 * t98 + t95 * t99) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t97, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:56
	% EndTime: 2019-10-09 22:21:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t128 = sin(pkin(6));
	t132 = sin(qJ(2));
	t141 = t128 * t132;
	t133 = cos(qJ(3));
	t140 = t128 * t133;
	t130 = cos(pkin(6));
	t139 = t130 * t132;
	t134 = cos(qJ(2));
	t138 = t130 * t134;
	t137 = qJD(2) * t133;
	t127 = sin(pkin(10));
	t129 = cos(pkin(10));
	t136 = t127 * t134 + t129 * t139;
	t135 = -t127 * t139 + t129 * t134;
	t131 = sin(qJ(3));
	t1 = [0, 0, t135 * qJD(2), 0, (t127 * t140 - t135 * t131) * qJD(3) + (-t127 * t138 - t129 * t132) * t137, 0; 0, 0, t136 * qJD(2), 0, (-t129 * t140 - t136 * t131) * qJD(3) + (-t127 * t132 + t129 * t138) * t137, 0; 0, 0, qJD(2) * t141, 0, t128 * t134 * t137 + (t130 * t133 - t131 * t141) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:56
	% EndTime: 2019-10-09 22:21:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t137 = sin(pkin(6));
	t141 = sin(qJ(2));
	t150 = t137 * t141;
	t142 = cos(qJ(3));
	t149 = t137 * t142;
	t139 = cos(pkin(6));
	t148 = t139 * t141;
	t143 = cos(qJ(2));
	t147 = t139 * t143;
	t146 = qJD(2) * t142;
	t136 = sin(pkin(10));
	t138 = cos(pkin(10));
	t145 = t136 * t143 + t138 * t148;
	t144 = -t136 * t148 + t138 * t143;
	t140 = sin(qJ(3));
	t1 = [0, 0, t144 * qJD(2), 0, (t136 * t149 - t144 * t140) * qJD(3) + (-t136 * t147 - t138 * t141) * t146, 0; 0, 0, t145 * qJD(2), 0, (-t138 * t149 - t145 * t140) * qJD(3) + (-t136 * t141 + t138 * t147) * t146, 0; 0, 0, qJD(2) * t150, 0, t137 * t143 * t146 + (t139 * t142 - t140 * t150) * qJD(3), 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end