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
% Datum: 2019-10-24 10:45
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:54
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
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(9);
	t13 = qJD(1) * cos(t12);
	t10 = qJD(1) * sin(t12);
	t1 = [0, 0, 0, 0, 0; t10, 0, 0, 0, 0; -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t13, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:53
	% EndTime: 2019-10-24 10:45:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t27 = qJD(1) + qJD(3);
	t26 = qJ(1) + pkin(9) + qJ(3);
	t24 = t27 * cos(t26);
	t23 = t27 * sin(t26);
	t1 = [0, 0, 0, 0, 0; t23, 0, t23, 0, 0; -t24, 0, -t24, 0, 0; 0, 0, 0, 0, 0; t24, 0, t24, 0, 0; t23, 0, t23, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:54
	% EndTime: 2019-10-24 10:45:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (90->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t116 = qJ(1) + pkin(9) + qJ(3);
	t114 = sin(t116);
	t117 = qJD(1) + qJD(3);
	t125 = t117 * t114;
	t115 = cos(t116);
	t124 = t117 * t115;
	t118 = sin(qJ(4));
	t123 = t117 * t118;
	t119 = cos(qJ(4));
	t122 = t117 * t119;
	t121 = qJD(4) * t118;
	t120 = qJD(4) * t119;
	t111 = -t114 * t121 + t115 * t122;
	t110 = t114 * t120 + t115 * t123;
	t109 = t114 * t122 + t115 * t121;
	t108 = t114 * t123 - t115 * t120;
	t1 = [0, 0, 0, -t121, 0; t109, 0, t109, t110, 0; -t111, 0, -t111, t108, 0; 0, 0, 0, -t120, 0; -t108, 0, -t108, t111, 0; t110, 0, t110, t109, 0; 0, 0, 0, 0, 0; -t124, 0, -t124, 0, 0; -t125, 0, -t125, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:45:54
	% EndTime: 2019-10-24 10:45:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (172->20), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t167 = qJ(4) + qJ(5);
	t163 = sin(t167);
	t165 = qJD(4) + qJD(5);
	t171 = t165 * t163;
	t164 = cos(t167);
	t170 = t165 * t164;
	t162 = qJ(1) + pkin(9) + qJ(3);
	t160 = sin(t162);
	t166 = qJD(1) + qJD(3);
	t169 = t166 * t160;
	t161 = cos(t162);
	t168 = t166 * t161;
	t157 = -t160 * t171 + t164 * t168;
	t156 = t160 * t170 + t163 * t168;
	t155 = t161 * t171 + t164 * t169;
	t154 = -t161 * t170 + t163 * t169;
	t1 = [0, 0, 0, -t171, -t171; t155, 0, t155, t156, t156; -t157, 0, -t157, t154, t154; 0, 0, 0, -t170, -t170; -t154, 0, -t154, t157, t157; t156, 0, t156, t155, t155; 0, 0, 0, 0, 0; -t168, 0, -t168, 0, 0; -t169, 0, -t169, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end