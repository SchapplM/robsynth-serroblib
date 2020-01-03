% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
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
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
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
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
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
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (138->10), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t130 = qJD(1) + qJD(2) + qJD(3);
	t132 = sin(qJ(4));
	t137 = t130 * t132;
	t133 = cos(qJ(4));
	t136 = t130 * t133;
	t135 = qJD(4) * t132;
	t134 = qJD(4) * t133;
	t131 = qJ(1) + qJ(2) + qJ(3);
	t129 = cos(t131);
	t128 = sin(t131);
	t127 = t130 * t129;
	t126 = t130 * t128;
	t125 = -t128 * t135 + t129 * t136;
	t124 = -t128 * t134 - t129 * t137;
	t123 = -t128 * t136 - t129 * t135;
	t122 = t128 * t137 - t129 * t134;
	t1 = [0, 0, 0, -t135, 0; t123, t123, t123, t124, 0; t125, t125, t125, -t122, 0; 0, 0, 0, -t134, 0; t122, t122, t122, -t125, 0; t124, t124, t124, t123, 0; 0, 0, 0, 0, 0; t127, t127, t127, 0, 0; t126, t126, t126, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (236->16), mult. (90->14), div. (0->0), fcn. (90->4), ass. (0->19)
	t162 = qJD(1) + qJD(2) + qJD(3);
	t167 = qJ(4) + qJ(5);
	t163 = sin(t167);
	t171 = t162 * t163;
	t164 = cos(t167);
	t170 = t162 * t164;
	t166 = qJD(4) + qJD(5);
	t169 = t166 * t163;
	t168 = t166 * t164;
	t165 = qJ(1) + qJ(2) + qJ(3);
	t161 = cos(t165);
	t160 = sin(t165);
	t159 = t162 * t161;
	t158 = t162 * t160;
	t155 = -t160 * t169 + t161 * t170;
	t154 = -t160 * t168 - t161 * t171;
	t153 = -t160 * t170 - t161 * t169;
	t152 = t160 * t171 - t161 * t168;
	t1 = [0, 0, 0, -t169, -t169; t153, t153, t153, t154, t154; t155, t155, t155, -t152, -t152; 0, 0, 0, -t168, -t168; t152, t152, t152, -t155, -t155; t154, t154, t154, t153, t153; 0, 0, 0, 0, 0; t159, t159, t159, 0, 0; t158, t158, t158, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end