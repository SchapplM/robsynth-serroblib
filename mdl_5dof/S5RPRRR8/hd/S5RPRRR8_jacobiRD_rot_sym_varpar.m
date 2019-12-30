% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:56
	% EndTime: 2019-12-29 17:48:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t14 = qJD(1) * sin(qJ(1));
	t13 = qJD(1) * cos(qJ(1));
	t1 = [-t13, 0, 0, 0, 0; -t14, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:56
	% EndTime: 2019-12-29 17:48:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->10), mult. (64->12), div. (0->0), fcn. (64->4), ass. (0->9)
	t95 = sin(qJ(1));
	t99 = qJD(3) * t95;
	t97 = cos(qJ(1));
	t98 = qJD(3) * t97;
	t96 = cos(qJ(3));
	t94 = sin(qJ(3));
	t89 = -t94 * t98 + t96 * t99 + (t94 * t97 - t95 * t96) * qJD(1);
	t88 = -t94 * t99 - t96 * t98 + (t94 * t95 + t96 * t97) * qJD(1);
	t1 = [-t88, 0, t88, 0, 0; t89, 0, -t89, 0, 0; 0, 0, 0, 0, 0; t89, 0, -t89, 0, 0; t88, 0, -t88, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:49:01
	% EndTime: 2019-12-29 17:49:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (72->15), mult. (190->16), div. (0->0), fcn. (202->6), ass. (0->18)
	t108 = qJD(1) - qJD(3);
	t104 = sin(qJ(3));
	t105 = sin(qJ(1));
	t106 = cos(qJ(3));
	t107 = cos(qJ(1));
	t88 = t107 * t104 - t105 * t106;
	t87 = -t105 * t104 - t107 * t106;
	t96 = sin(qJ(4));
	t103 = qJD(4) * t96;
	t97 = cos(qJ(4));
	t102 = qJD(4) * t97;
	t85 = t108 * t87;
	t81 = t88 * t102 + t85 * t96;
	t82 = t88 * t103 - t85 * t97;
	t86 = t108 * t88;
	t83 = t87 * t102 - t86 * t96;
	t84 = t87 * t103 + t86 * t97;
	t1 = [-t82, 0, t82, t83, 0; t84, 0, -t84, t81, 0; 0, 0, 0, t103, 0; -t81, 0, t81, -t84, 0; t83, 0, -t83, -t82, 0; 0, 0, 0, t102, 0; -t86, 0, t86, 0, 0; t85, 0, -t85, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:56
	% EndTime: 2019-12-29 17:48:56
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (148->22), mult. (244->22), div. (0->0), fcn. (260->6), ass. (0->21)
	t133 = sin(qJ(3));
	t134 = sin(qJ(1));
	t135 = cos(qJ(3));
	t136 = cos(qJ(1));
	t122 = t136 * t133 - t134 * t135;
	t132 = qJ(4) + qJ(5);
	t129 = sin(t132);
	t131 = qJD(4) + qJD(5);
	t123 = t131 * t129;
	t130 = cos(t132);
	t124 = t131 * t130;
	t139 = qJD(3) * t134;
	t138 = qJD(3) * t136;
	t137 = t134 * t133 + t136 * t135;
	t117 = t137 * qJD(1) - t133 * t139 - t135 * t138;
	t113 = -t117 * t129 + t122 * t124;
	t114 = t117 * t130 + t122 * t123;
	t118 = t122 * qJD(1) - t133 * t138 + t135 * t139;
	t115 = -t118 * t129 - t124 * t137;
	t116 = t118 * t130 - t123 * t137;
	t1 = [-t114, 0, t114, t115, t115; t116, 0, -t116, t113, t113; 0, 0, 0, t123, t123; -t113, 0, t113, -t116, -t116; t115, 0, -t115, -t114, -t114; 0, 0, 0, t124, t124; -t118, 0, t118, 0, 0; -t117, 0, t117, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end