% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:49
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:06
	% EndTime: 2019-10-24 10:49:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:06
	% EndTime: 2019-10-24 10:49:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:06
	% EndTime: 2019-10-24 10:49:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t26 = qJ(1) + qJ(2);
	t25 = qJD(1) + qJD(2);
	t23 = t25 * cos(t26);
	t22 = t25 * sin(t26);
	t1 = [0, 0, 0, 0, 0; t22, t22, 0, 0, 0; -t23, -t23, 0, 0, 0; 0, 0, 0, 0, 0; t23, t23, 0, 0, 0; t22, t22, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:07
	% EndTime: 2019-10-24 10:49:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t31 = qJD(1) + qJD(2);
	t30 = qJ(1) + qJ(2) + pkin(9);
	t28 = t31 * cos(t30);
	t27 = t31 * sin(t30);
	t1 = [0, 0, 0, 0, 0; t27, t27, 0, 0, 0; -t28, -t28, 0, 0, 0; 0, 0, 0, 0, 0; t28, t28, 0, 0, 0; t27, t27, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:07
	% EndTime: 2019-10-24 10:49:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (90->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t121 = qJ(1) + qJ(2) + pkin(9);
	t119 = sin(t121);
	t122 = qJD(1) + qJD(2);
	t130 = t122 * t119;
	t120 = cos(t121);
	t129 = t122 * t120;
	t123 = sin(qJ(4));
	t128 = t122 * t123;
	t124 = cos(qJ(4));
	t127 = t122 * t124;
	t126 = qJD(4) * t123;
	t125 = qJD(4) * t124;
	t116 = -t119 * t126 + t120 * t127;
	t115 = t119 * t125 + t120 * t128;
	t114 = t119 * t127 + t120 * t126;
	t113 = t119 * t128 - t120 * t125;
	t1 = [0, 0, 0, -t126, 0; t114, t114, 0, t115, 0; -t116, -t116, 0, t113, 0; 0, 0, 0, -t125, 0; -t113, -t113, 0, t116, 0; t115, t115, 0, t114, 0; 0, 0, 0, 0, 0; -t129, -t129, 0, 0, 0; -t130, -t130, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:07
	% EndTime: 2019-10-24 10:49:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (172->20), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t170 = qJ(4) + qJ(5);
	t166 = sin(t170);
	t168 = qJD(4) + qJD(5);
	t174 = t168 * t166;
	t167 = cos(t170);
	t173 = t168 * t167;
	t165 = qJ(1) + qJ(2) + pkin(9);
	t163 = sin(t165);
	t169 = qJD(1) + qJD(2);
	t172 = t169 * t163;
	t164 = cos(t165);
	t171 = t169 * t164;
	t160 = -t163 * t174 + t167 * t171;
	t159 = t163 * t173 + t166 * t171;
	t158 = t164 * t174 + t167 * t172;
	t157 = -t164 * t173 + t166 * t172;
	t1 = [0, 0, 0, -t174, -t174; t158, t158, 0, t159, t159; -t160, -t160, 0, t157, t157; 0, 0, 0, -t173, -t173; -t157, -t157, 0, t160, t160; t159, t159, 0, t158, t158; 0, 0, 0, 0, 0; -t171, -t171, 0, 0, 0; -t172, -t172, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end