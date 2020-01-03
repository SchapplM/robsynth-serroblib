% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
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
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:12
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
	% StartTime: 2020-01-03 11:40:12
	% EndTime: 2020-01-03 11:40:12
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
	% StartTime: 2020-01-03 11:40:12
	% EndTime: 2020-01-03 11:40:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->10), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t102 = qJ(1) + pkin(8);
	t98 = sin(t102);
	t107 = qJD(1) * t98;
	t101 = qJ(3) + pkin(9);
	t97 = sin(t101);
	t106 = qJD(3) * t97;
	t99 = cos(t101);
	t105 = qJD(3) * t99;
	t100 = cos(t102);
	t104 = qJD(1) * t100;
	t103 = qJD(3) * t100;
	t96 = t99 * t104 - t98 * t106;
	t95 = -t97 * t104 - t98 * t105;
	t94 = -t97 * t103 - t99 * t107;
	t93 = -t99 * t103 + t97 * t107;
	t1 = [0, 0, -t106, 0, 0; t94, 0, t95, 0, 0; t96, 0, -t93, 0, 0; 0, 0, -t105, 0, 0; t93, 0, -t96, 0, 0; t95, 0, t94, 0, 0; 0, 0, 0, 0, 0; t104, 0, 0, 0, 0; t107, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:12
	% EndTime: 2020-01-03 11:40:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (114->15), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t137 = qJ(3) + pkin(9) + qJ(5);
	t133 = sin(t137);
	t138 = qJD(3) + qJD(5);
	t143 = t138 * t133;
	t134 = cos(t137);
	t142 = t138 * t134;
	t139 = qJ(1) + pkin(8);
	t135 = sin(t139);
	t141 = qJD(1) * t135;
	t136 = cos(t139);
	t140 = qJD(1) * t136;
	t130 = t134 * t140 - t135 * t143;
	t129 = -t133 * t140 - t135 * t142;
	t128 = -t134 * t141 - t136 * t143;
	t127 = t133 * t141 - t136 * t142;
	t1 = [0, 0, -t143, 0, -t143; t128, 0, t129, 0, t129; t130, 0, -t127, 0, -t127; 0, 0, -t142, 0, -t142; t127, 0, -t130, 0, -t130; t129, 0, t128, 0, t128; 0, 0, 0, 0, 0; t140, 0, 0, 0, 0; t141, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end