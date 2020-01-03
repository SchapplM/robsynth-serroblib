% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:14
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
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:15
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
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t49 = qJD(1) * sin(pkin(9));
	t48 = qJD(1) * cos(pkin(9));
	t45 = qJ(1) + pkin(8);
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [0, 0, 0, 0, 0; -t43 * t48, 0, 0, 0, 0; t44 * t48, 0, 0, 0, 0; 0, 0, 0, 0, 0; t43 * t49, 0, 0, 0, 0; -t44 * t49, 0, 0, 0, 0; 0, 0, 0, 0, 0; qJD(1) * t44, 0, 0, 0, 0; qJD(1) * t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:15
	% EndTime: 2020-01-03 11:29:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->10), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t95 = qJ(1) + pkin(8);
	t91 = sin(t95);
	t100 = qJD(1) * t91;
	t93 = cos(t95);
	t99 = qJD(1) * t93;
	t94 = pkin(9) + qJ(4);
	t90 = sin(t94);
	t98 = qJD(4) * t90;
	t92 = cos(t94);
	t97 = qJD(4) * t92;
	t96 = qJD(4) * t93;
	t89 = -t91 * t98 + t92 * t99;
	t88 = -t90 * t99 - t91 * t97;
	t87 = -t92 * t100 - t90 * t96;
	t86 = t90 * t100 - t92 * t96;
	t1 = [0, 0, 0, -t98, 0; t87, 0, 0, t88, 0; t89, 0, 0, -t86, 0; 0, 0, 0, -t97, 0; t86, 0, 0, -t89, 0; t88, 0, 0, t87, 0; 0, 0, 0, 0, 0; t99, 0, 0, 0, 0; t100, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:15
	% EndTime: 2020-01-03 11:29:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (114->15), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t132 = pkin(9) + qJ(4) + qJ(5);
	t128 = sin(t132);
	t133 = qJD(4) + qJD(5);
	t138 = t133 * t128;
	t129 = cos(t132);
	t137 = t133 * t129;
	t134 = qJ(1) + pkin(8);
	t130 = sin(t134);
	t136 = qJD(1) * t130;
	t131 = cos(t134);
	t135 = qJD(1) * t131;
	t125 = t129 * t135 - t130 * t138;
	t124 = -t128 * t135 - t130 * t137;
	t123 = -t129 * t136 - t131 * t138;
	t122 = t128 * t136 - t131 * t137;
	t1 = [0, 0, 0, -t138, -t138; t123, 0, 0, t124, t124; t125, 0, 0, -t122, -t122; 0, 0, 0, -t137, -t137; t122, 0, 0, -t125, -t125; t124, 0, 0, t123, t123; 0, 0, 0, 0, 0; t135, 0, 0, 0, 0; t136, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end