% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
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
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
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
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t118 = qJD(1) + qJD(2);
	t120 = sin(qJ(3));
	t125 = t118 * t120;
	t121 = cos(qJ(3));
	t124 = t118 * t121;
	t123 = qJD(3) * t120;
	t122 = qJD(3) * t121;
	t119 = qJ(1) + qJ(2);
	t117 = cos(t119);
	t116 = sin(t119);
	t115 = t118 * t117;
	t114 = t118 * t116;
	t113 = -t116 * t123 + t117 * t124;
	t112 = -t116 * t122 - t117 * t125;
	t111 = -t116 * t124 - t117 * t123;
	t110 = t116 * t125 - t117 * t122;
	t1 = [0, 0, -t123, 0, 0; t111, t111, t112, 0, 0; t113, t113, -t110, 0, 0; 0, 0, -t122, 0, 0; t110, t110, -t113, 0, 0; t112, t112, t111, 0, 0; 0, 0, 0, 0, 0; t115, t115, 0, 0, 0; t114, t114, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (132->16), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t152 = qJ(3) + qJ(4);
	t146 = sin(t152);
	t150 = qJD(3) + qJD(4);
	t155 = t150 * t146;
	t148 = cos(t152);
	t154 = t150 * t148;
	t153 = qJ(1) + qJ(2);
	t147 = sin(t153);
	t151 = qJD(1) + qJD(2);
	t144 = t151 * t147;
	t149 = cos(t153);
	t145 = t151 * t149;
	t141 = t148 * t145 - t147 * t155;
	t140 = -t146 * t145 - t147 * t154;
	t139 = -t148 * t144 - t149 * t155;
	t138 = t146 * t144 - t149 * t154;
	t1 = [0, 0, -t155, -t155, 0; t139, t139, t140, t140, 0; t141, t141, -t138, -t138, 0; 0, 0, -t154, -t154, 0; t138, t138, -t141, -t141, 0; t140, t140, t139, t139, 0; 0, 0, 0, 0, 0; t145, t145, 0, 0, 0; t144, t144, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (132->16), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t164 = qJ(3) + qJ(4);
	t158 = sin(t164);
	t162 = qJD(3) + qJD(4);
	t167 = t162 * t158;
	t160 = cos(t164);
	t166 = t162 * t160;
	t165 = qJ(1) + qJ(2);
	t159 = sin(t165);
	t163 = qJD(1) + qJD(2);
	t156 = t163 * t159;
	t161 = cos(t165);
	t157 = t163 * t161;
	t153 = t160 * t157 - t159 * t167;
	t152 = -t158 * t157 - t159 * t166;
	t151 = -t160 * t156 - t161 * t167;
	t150 = t158 * t156 - t161 * t166;
	t1 = [0, 0, -t167, -t167, 0; t151, t151, t152, t152, 0; t153, t153, -t150, -t150, 0; 0, 0, -t166, -t166, 0; t150, t150, -t153, -t153, 0; t152, t152, t151, t151, 0; 0, 0, 0, 0, 0; t157, t157, 0, 0, 0; t156, t156, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end