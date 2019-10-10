% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t78 = sin(qJ(1));
	t83 = qJD(1) * t78;
	t80 = cos(qJ(1));
	t82 = qJD(1) * t80;
	t81 = qJD(2) * cos(qJ(2));
	t77 = sin(qJ(2));
	t1 = [0, t82, 0, -t77 * t83 + t80 * t81, 0, 0; 0, t83, 0, t77 * t82 + t78 * t81, 0, 0; 0, 0, 0, qJD(2) * t77, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t93 = sin(qJ(1));
	t98 = qJD(1) * t93;
	t95 = cos(qJ(1));
	t97 = qJD(1) * t95;
	t96 = qJD(2) * cos(qJ(2));
	t92 = sin(qJ(2));
	t1 = [0, t97, 0, -t92 * t98 + t95 * t96, 0, 0; 0, t98, 0, t92 * t97 + t93 * t96, 0, 0; 0, 0, 0, qJD(2) * t92, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->10)
	t121 = sin(qJ(1));
	t127 = qJD(1) * t121;
	t123 = cos(qJ(1));
	t126 = qJD(1) * t123;
	t120 = sin(qJ(2));
	t125 = qJD(2) * t120;
	t124 = qJD(2) * cos(qJ(2));
	t119 = t120 * t126 + t121 * t124;
	t118 = t120 * t127 - t123 * t124;
	t1 = [0, t126, 0, -t118, 0, t118; 0, t127, 0, t119, 0, -t119; 0, 0, 0, t125, 0, -t125;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end