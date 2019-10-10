% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:58
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t30 = qJD(1) * cos(qJ(1));
	t29 = qJD(1) * sin(qJ(1));
	t1 = [0, t30, t30, 0, 0, 0; 0, t29, t29, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (6->2), div. (0->0), fcn. (6->2), ass. (0->3)
	t36 = qJD(1) * cos(qJ(1));
	t35 = qJD(1) * sin(qJ(1));
	t1 = [0, t36, t36, t36, 0, 0; 0, t35, t35, t35, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->4), mult. (15->8), div. (0->0), fcn. (15->4), ass. (0->9)
	t120 = qJD(2) + qJD(3) + qJD(4);
	t123 = qJ(2) + qJ(3) + qJ(4);
	t126 = cos(t123) * t120;
	t124 = sin(qJ(1));
	t121 = qJD(1) * t124;
	t125 = cos(qJ(1));
	t122 = qJD(1) * t125;
	t118 = sin(t123);
	t1 = [0, t122, t122, t122, -t118 * t121 + t125 * t126, 0; 0, t121, t121, t121, t118 * t122 + t124 * t126, 0; 0, 0, 0, 0, t120 * t118, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->4), mult. (15->8), div. (0->0), fcn. (15->4), ass. (0->9)
	t135 = qJD(2) + qJD(3) + qJD(4);
	t138 = qJ(2) + qJ(3) + qJ(4);
	t141 = cos(t138) * t135;
	t139 = sin(qJ(1));
	t136 = qJD(1) * t139;
	t140 = cos(qJ(1));
	t137 = qJD(1) * t140;
	t133 = sin(t138);
	t1 = [0, t137, t137, t137, -t133 * t136 + t140 * t141, 0; 0, t136, t136, t136, t133 * t137 + t139 * t141, 0; 0, 0, 0, 0, t135 * t133, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end