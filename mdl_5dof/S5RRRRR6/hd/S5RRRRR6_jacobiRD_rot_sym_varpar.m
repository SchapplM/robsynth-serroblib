% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR6
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
% Datum: 2019-10-24 10:52
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
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
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:49
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
	% StartTime: 2019-10-24 10:52:49
	% EndTime: 2019-10-24 10:52:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (62->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t116 = qJ(1) + qJ(2);
	t113 = sin(t116);
	t115 = qJD(1) + qJD(2);
	t124 = t115 * t113;
	t114 = cos(t116);
	t123 = t115 * t114;
	t117 = sin(qJ(3));
	t122 = t115 * t117;
	t118 = cos(qJ(3));
	t121 = t115 * t118;
	t120 = qJD(3) * t117;
	t119 = qJD(3) * t118;
	t110 = -t113 * t120 + t114 * t121;
	t109 = t113 * t119 + t114 * t122;
	t108 = t113 * t121 + t114 * t120;
	t107 = t113 * t122 - t114 * t119;
	t1 = [0, 0, -t120, 0, 0; t108, t108, t109, 0, 0; -t110, -t110, t107, 0, 0; 0, 0, -t119, 0, 0; -t107, -t107, t110, 0, 0; t109, t109, t108, 0, 0; 0, 0, 0, 0, 0; -t123, -t123, 0, 0, 0; -t124, -t124, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:50
	% EndTime: 2019-10-24 10:52:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (136->20), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t165 = qJ(3) + qJ(4);
	t159 = sin(t165);
	t163 = qJD(3) + qJD(4);
	t170 = t163 * t159;
	t161 = cos(t165);
	t169 = t163 * t161;
	t166 = qJ(1) + qJ(2);
	t160 = sin(t166);
	t164 = qJD(1) + qJD(2);
	t168 = t164 * t160;
	t162 = cos(t166);
	t167 = t164 * t162;
	t156 = -t160 * t170 + t161 * t167;
	t155 = t159 * t167 + t160 * t169;
	t154 = t161 * t168 + t162 * t170;
	t153 = t159 * t168 - t162 * t169;
	t1 = [0, 0, -t170, -t170, 0; t154, t154, t155, t155, 0; -t156, -t156, t153, t153, 0; 0, 0, -t169, -t169, 0; -t153, -t153, t156, t156, 0; t155, t155, t154, t154, 0; 0, 0, 0, 0, 0; -t167, -t167, 0, 0, 0; -t168, -t168, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:50
	% EndTime: 2019-10-24 10:52:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (242->22), mult. (90->14), div. (0->0), fcn. (90->4), ass. (0->19)
	t179 = qJ(3) + qJ(4) + qJ(5);
	t174 = sin(t179);
	t176 = qJD(3) + qJD(4) + qJD(5);
	t187 = t176 * t174;
	t175 = cos(t179);
	t186 = t176 * t175;
	t181 = qJ(1) + qJ(2);
	t177 = sin(t181);
	t185 = t176 * t177;
	t178 = cos(t181);
	t184 = t176 * t178;
	t180 = qJD(1) + qJD(2);
	t183 = t180 * t177;
	t182 = t180 * t178;
	t171 = -t174 * t185 + t175 * t182;
	t170 = t174 * t182 + t175 * t185;
	t169 = t174 * t184 + t175 * t183;
	t168 = t174 * t183 - t175 * t184;
	t1 = [0, 0, -t187, -t187, -t187; t169, t169, t170, t170, t170; -t171, -t171, t168, t168, t168; 0, 0, -t186, -t186, -t186; -t168, -t168, t171, t171, t171; t170, t170, t169, t169, t169; 0, 0, 0, 0, 0; -t182, -t182, 0, 0, 0; -t183, -t183, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end