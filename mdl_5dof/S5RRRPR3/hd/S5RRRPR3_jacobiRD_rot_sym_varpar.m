% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR3
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
% Datum: 2019-10-24 10:50
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
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
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
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
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
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
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (88->17), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t147 = qJ(1) + qJ(2);
	t143 = sin(t147);
	t145 = qJD(1) + qJD(2);
	t151 = t145 * t143;
	t144 = cos(t147);
	t150 = t145 * t144;
	t149 = qJD(3) * t143;
	t148 = qJD(3) * t144;
	t146 = qJ(3) + pkin(9);
	t142 = cos(t146);
	t141 = sin(t146);
	t138 = -t141 * t149 + t142 * t150;
	t137 = t141 * t150 + t142 * t149;
	t136 = t141 * t148 + t142 * t151;
	t135 = t141 * t151 - t142 * t148;
	t1 = [0, 0, -qJD(3) * t141, 0, 0; t136, t136, t137, 0, 0; -t138, -t138, t135, 0, 0; 0, 0, -qJD(3) * t142, 0, 0; -t135, -t135, t138, 0, 0; t137, t137, t136, 0, 0; 0, 0, 0, 0, 0; -t150, -t150, 0, 0, 0; -t151, -t151, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (172->20), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t166 = qJ(3) + pkin(9) + qJ(5);
	t164 = sin(t166);
	t169 = qJD(3) + qJD(5);
	t175 = t169 * t164;
	t165 = cos(t166);
	t174 = t169 * t165;
	t171 = qJ(1) + qJ(2);
	t167 = sin(t171);
	t170 = qJD(1) + qJD(2);
	t173 = t170 * t167;
	t168 = cos(t171);
	t172 = t170 * t168;
	t161 = t165 * t172 - t167 * t175;
	t160 = t164 * t172 + t167 * t174;
	t159 = t165 * t173 + t168 * t175;
	t158 = t164 * t173 - t168 * t174;
	t1 = [0, 0, -t175, 0, -t175; t159, t159, t160, 0, t160; -t161, -t161, t158, 0, t158; 0, 0, -t174, 0, -t174; -t158, -t158, t161, 0, t161; t160, t160, t159, 0, t159; 0, 0, 0, 0, 0; -t172, -t172, 0, 0, 0; -t173, -t173, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end