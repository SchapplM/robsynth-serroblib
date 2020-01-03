% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR5
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(1) * sin(qJ(1));
	t9 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
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
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.05s
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
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (86->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t120 = qJD(1) + qJD(3);
	t121 = sin(qJ(4));
	t126 = t120 * t121;
	t122 = cos(qJ(4));
	t125 = t120 * t122;
	t124 = qJD(4) * t121;
	t123 = qJD(4) * t122;
	t119 = qJ(1) + pkin(9) + qJ(3);
	t118 = cos(t119);
	t117 = sin(t119);
	t116 = t120 * t118;
	t115 = t120 * t117;
	t114 = -t117 * t124 + t118 * t125;
	t113 = -t117 * t123 - t118 * t126;
	t112 = -t117 * t125 - t118 * t124;
	t111 = t117 * t126 - t118 * t123;
	t1 = [0, 0, 0, -t124, 0; t112, 0, t112, t113, 0; t114, 0, t114, -t111, 0; 0, 0, 0, -t123, 0; t111, 0, t111, -t114, 0; t113, 0, t113, t112, 0; 0, 0, 0, 0, 0; t116, 0, t116, 0, 0; t115, 0, t115, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (168->16), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t156 = qJ(4) + qJ(5);
	t152 = sin(t156);
	t154 = qJD(4) + qJD(5);
	t158 = t154 * t152;
	t153 = cos(t156);
	t157 = t154 * t153;
	t151 = qJ(1) + pkin(9) + qJ(3);
	t149 = sin(t151);
	t155 = qJD(1) + qJD(3);
	t147 = t155 * t149;
	t150 = cos(t151);
	t148 = t155 * t150;
	t144 = t153 * t148 - t149 * t158;
	t143 = -t152 * t148 - t149 * t157;
	t142 = -t153 * t147 - t150 * t158;
	t141 = t152 * t147 - t150 * t157;
	t1 = [0, 0, 0, -t158, -t158; t142, 0, t142, t143, t143; t144, 0, t144, -t141, -t141; 0, 0, 0, -t157, -t157; t141, 0, t141, -t144, -t144; t143, 0, t143, t142, t142; 0, 0, 0, 0, 0; t148, 0, t148, 0, 0; t147, 0, t147, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end