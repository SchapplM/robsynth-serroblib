% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
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
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t23 = qJD(1) + qJD(2);
	t24 = qJ(1) + qJ(2);
	t25 = t23 * sin(t24);
	t20 = t23 * cos(t24);
	t1 = [0, 0, 0, 0, 0; -t25, -t25, 0, 0, 0; t20, t20, 0, 0, 0; 0, 0, 0, 0, 0; -t20, -t20, 0, 0, 0; -t25, -t25, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (57->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t35 = qJD(1) + qJD(2) + qJD(3);
	t36 = qJ(1) + qJ(2) + qJ(3);
	t37 = t35 * sin(t36);
	t32 = t35 * cos(t36);
	t1 = [0, 0, 0, 0, 0; -t37, -t37, -t37, 0, 0; t32, t32, t32, 0, 0; 0, 0, 0, 0, 0; -t32, -t32, -t32, 0, 0; -t37, -t37, -t37, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (69->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t41 = qJ(1) + qJ(2) + qJ(3) + pkin(9);
	t42 = qJD(1) + qJD(2) + qJD(3);
	t43 = t42 * sin(t41);
	t38 = t42 * cos(t41);
	t1 = [0, 0, 0, 0, 0; -t43, -t43, -t43, 0, 0; t38, t38, t38, 0, 0; 0, 0, 0, 0, 0; -t38, -t38, -t38, 0, 0; -t43, -t43, -t43, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (176->10), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t135 = qJD(1) + qJD(2) + qJD(3);
	t136 = sin(qJ(5));
	t141 = t135 * t136;
	t137 = cos(qJ(5));
	t140 = t135 * t137;
	t139 = qJD(5) * t136;
	t138 = qJD(5) * t137;
	t134 = qJ(1) + qJ(2) + qJ(3) + pkin(9);
	t133 = cos(t134);
	t132 = sin(t134);
	t131 = t135 * t133;
	t130 = t135 * t132;
	t129 = -t132 * t139 + t133 * t140;
	t128 = -t132 * t138 - t133 * t141;
	t127 = -t132 * t140 - t133 * t139;
	t126 = t132 * t141 - t133 * t138;
	t1 = [0, 0, 0, 0, -t139; t127, t127, t127, 0, t128; t129, t129, t129, 0, -t126; 0, 0, 0, 0, -t138; t126, t126, t126, 0, -t129; t128, t128, t128, 0, t127; 0, 0, 0, 0, 0; t131, t131, t131, 0, 0; t130, t130, t130, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end