% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:50
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:51
	% EndTime: 2020-01-03 11:52:51
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
	% StartTime: 2020-01-03 11:52:50
	% EndTime: 2020-01-03 11:52:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(9);
	t14 = qJD(1) * sin(t12);
	t13 = qJD(1) * cos(t12);
	t1 = [0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; -t14, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:51
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t24 = qJ(1) + pkin(9) + qJ(3);
	t25 = qJD(1) + qJD(3);
	t26 = t25 * sin(t24);
	t21 = t25 * cos(t24);
	t1 = [0, 0, 0, 0, 0; -t26, 0, -t26, 0, 0; t21, 0, t21, 0, 0; 0, 0, 0, 0, 0; -t21, 0, -t21, 0, 0; -t26, 0, -t26, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:51
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (69->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t36 = qJ(1) + pkin(9) + qJ(3) + qJ(4);
	t37 = qJD(1) + qJD(3) + qJD(4);
	t38 = t37 * sin(t36);
	t33 = t37 * cos(t36);
	t1 = [0, 0, 0, 0, 0; -t38, 0, -t38, -t38, 0; t33, 0, t33, t33, 0; 0, 0, 0, 0, 0; -t33, 0, -t33, -t33, 0; -t38, 0, -t38, -t38, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:51
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (176->10), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t132 = qJD(1) + qJD(3) + qJD(4);
	t133 = sin(qJ(5));
	t138 = t132 * t133;
	t134 = cos(qJ(5));
	t137 = t132 * t134;
	t136 = qJD(5) * t133;
	t135 = qJD(5) * t134;
	t131 = qJ(1) + pkin(9) + qJ(3) + qJ(4);
	t130 = cos(t131);
	t129 = sin(t131);
	t128 = t132 * t130;
	t127 = t132 * t129;
	t126 = -t129 * t136 + t130 * t137;
	t125 = -t129 * t135 - t130 * t138;
	t124 = -t129 * t137 - t130 * t136;
	t123 = t129 * t138 - t130 * t135;
	t1 = [0, 0, 0, 0, -t136; t124, 0, t124, t124, t125; t126, 0, t126, t126, -t123; 0, 0, 0, 0, -t135; t123, 0, t123, t123, -t126; t125, 0, t125, t125, t124; 0, 0, 0, 0, 0; t128, 0, t128, t128, 0; t127, 0, t127, t127, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end