% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:40
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR7_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t98 = cos(qJ(2));
	t99 = cos(pkin(6)) * t98;
	t97 = sin(qJ(2));
	t95 = cos(pkin(10));
	t94 = sin(pkin(10));
	t1 = [0, 0, 0, (-t94 * t99 - t95 * t97) * qJD(2), 0, 0; 0, 0, 0, (-t94 * t97 + t95 * t99) * qJD(2), 0, 0; 0, 0, 0, sin(pkin(6)) * qJD(2) * t98, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (12->11), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t126 = sin(pkin(6));
	t129 = sin(qJ(4));
	t139 = t126 * t129;
	t132 = cos(qJ(2));
	t138 = t126 * t132;
	t128 = cos(pkin(6));
	t130 = sin(qJ(2));
	t137 = t128 * t130;
	t136 = t128 * t132;
	t135 = qJD(2) * t129;
	t125 = sin(pkin(10));
	t127 = cos(pkin(10));
	t134 = -t125 * t130 + t127 * t136;
	t133 = t125 * t136 + t127 * t130;
	t131 = cos(qJ(4));
	t1 = [0, 0, 0, -t133 * qJD(2), 0, (-t125 * t139 + t133 * t131) * qJD(4) + (-t125 * t137 + t127 * t132) * t135; 0, 0, 0, t134 * qJD(2), 0, (t127 * t139 - t134 * t131) * qJD(4) + (t125 * t132 + t127 * t137) * t135; 0, 0, 0, qJD(2) * t138, 0, t126 * t130 * t135 + (-t128 * t129 - t131 * t138) * qJD(4);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end