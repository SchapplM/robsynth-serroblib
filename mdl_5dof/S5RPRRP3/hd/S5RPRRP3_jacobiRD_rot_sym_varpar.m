% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(1) * sin(qJ(1));
	t9 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t14 = qJD(1) * sin(t12);
	t13 = qJD(1) * cos(t12);
	t1 = [0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; -t14, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->9), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t84 = sin(qJ(3));
	t89 = qJD(1) * t84;
	t85 = cos(qJ(3));
	t88 = qJD(1) * t85;
	t87 = qJD(3) * t84;
	t86 = qJD(3) * t85;
	t83 = qJ(1) + pkin(8);
	t82 = cos(t83);
	t81 = sin(t83);
	t80 = -t81 * t87 + t82 * t88;
	t79 = -t81 * t86 - t82 * t89;
	t78 = -t81 * t88 - t82 * t87;
	t77 = t81 * t89 - t82 * t86;
	t1 = [0, 0, -t87, 0, 0; t78, 0, t79, 0, 0; t80, 0, -t77, 0, 0; 0, 0, -t86, 0, 0; t77, 0, -t80, 0, 0; t79, 0, t78, 0, 0; 0, 0, 0, 0, 0; qJD(1) * t82, 0, 0, 0, 0; qJD(1) * t81, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (86->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t130 = qJ(3) + qJ(4);
	t126 = sin(t130);
	t128 = qJD(3) + qJD(4);
	t134 = t128 * t126;
	t127 = cos(t130);
	t133 = t128 * t127;
	t132 = qJD(1) * t126;
	t131 = qJD(1) * t127;
	t129 = qJ(1) + pkin(8);
	t125 = cos(t129);
	t124 = sin(t129);
	t121 = -t124 * t134 + t125 * t131;
	t120 = -t124 * t133 - t125 * t132;
	t119 = -t124 * t131 - t125 * t134;
	t118 = t124 * t132 - t125 * t133;
	t1 = [0, 0, -t134, -t134, 0; t119, 0, t120, t120, 0; t121, 0, -t118, -t118, 0; 0, 0, -t133, -t133, 0; t118, 0, -t121, -t121, 0; t120, 0, t119, t119, 0; 0, 0, 0, 0, 0; qJD(1) * t125, 0, 0, 0, 0; qJD(1) * t124, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (86->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t137 = qJ(3) + qJ(4);
	t133 = sin(t137);
	t135 = qJD(3) + qJD(4);
	t141 = t135 * t133;
	t134 = cos(t137);
	t140 = t135 * t134;
	t139 = qJD(1) * t133;
	t138 = qJD(1) * t134;
	t136 = qJ(1) + pkin(8);
	t132 = cos(t136);
	t131 = sin(t136);
	t128 = -t131 * t141 + t132 * t138;
	t127 = -t131 * t140 - t132 * t139;
	t126 = -t131 * t138 - t132 * t141;
	t125 = t131 * t139 - t132 * t140;
	t1 = [0, 0, -t141, -t141, 0; t126, 0, t127, t127, 0; t128, 0, -t125, -t125, 0; 0, 0, -t140, -t140, 0; t125, 0, -t128, -t128, 0; t127, 0, t126, t126, 0; 0, 0, 0, 0, 0; qJD(1) * t132, 0, 0, 0, 0; qJD(1) * t131, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end