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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 18:48:50
	% EndTime: 2019-12-05 18:48:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:48:50
	% EndTime: 2019-12-05 18:48:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:48:50
	% EndTime: 2019-12-05 18:48:50
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-12-05 18:48:50
	% EndTime: 2019-12-05 18:48:50
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
	% StartTime: 2019-12-05 18:48:50
	% EndTime: 2019-12-05 18:48:51
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
	% StartTime: 2019-12-05 18:48:51
	% EndTime: 2019-12-05 18:48:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (136->20), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t166 = qJ(3) + qJ(4);
	t160 = sin(t166);
	t164 = qJD(3) + qJD(4);
	t171 = t164 * t160;
	t162 = cos(t166);
	t170 = t164 * t162;
	t167 = qJ(1) + qJ(2);
	t161 = sin(t167);
	t165 = qJD(1) + qJD(2);
	t169 = t165 * t161;
	t163 = cos(t167);
	t168 = t165 * t163;
	t157 = -t161 * t171 + t162 * t168;
	t156 = t160 * t168 + t161 * t170;
	t155 = t162 * t169 + t163 * t171;
	t154 = t160 * t169 - t163 * t170;
	t1 = [0, 0, -t171, -t171, 0; t155, t155, t156, t156, 0; -t157, -t157, t154, t154, 0; 0, 0, -t170, -t170, 0; -t154, -t154, t157, t157, 0; t156, t156, t155, t155, 0; 0, 0, 0, 0, 0; -t168, -t168, 0, 0, 0; -t169, -t169, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end