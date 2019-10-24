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
% Datum: 2019-10-24 10:52
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:52:25
	% EndTime: 2019-10-24 10:52:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
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
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
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
	% StartTime: 2019-10-24 10:52:25
	% EndTime: 2019-10-24 10:52:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (51->5), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t38 = qJ(1) + qJ(2) + qJ(3);
	t37 = qJD(1) + qJD(2) + qJD(3);
	t35 = t37 * cos(t38);
	t34 = t37 * sin(t38);
	t1 = [0, 0, 0, 0, 0; t34, t34, t34, 0, 0; -t35, -t35, -t35, 0, 0; 0, 0, 0, 0, 0; t35, t35, t35, 0, 0; t34, t34, t34, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (144->20), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t128 = qJ(1) + qJ(2) + qJ(3);
	t125 = sin(t128);
	t127 = qJD(1) + qJD(2) + qJD(3);
	t136 = t127 * t125;
	t126 = cos(t128);
	t135 = t127 * t126;
	t129 = sin(qJ(4));
	t134 = t127 * t129;
	t130 = cos(qJ(4));
	t133 = t127 * t130;
	t132 = qJD(4) * t129;
	t131 = qJD(4) * t130;
	t122 = -t125 * t132 + t126 * t133;
	t121 = t125 * t131 + t126 * t134;
	t120 = t125 * t133 + t126 * t132;
	t119 = t125 * t134 - t126 * t131;
	t1 = [0, 0, 0, -t132, 0; t120, t120, t120, t121, 0; -t122, -t122, -t122, t119, 0; 0, 0, 0, -t131, 0; -t119, -t119, -t119, t122, 0; t121, t121, t121, t120, 0; 0, 0, 0, 0, 0; -t135, -t135, -t135, 0, 0; -t136, -t136, -t136, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (242->24), mult. (90->14), div. (0->0), fcn. (90->4), ass. (0->19)
	t176 = qJ(1) + qJ(2) + qJ(3);
	t171 = sin(t176);
	t173 = qJD(1) + qJD(2) + qJD(3);
	t184 = t173 * t171;
	t172 = cos(t176);
	t183 = t173 * t172;
	t178 = qJ(4) + qJ(5);
	t174 = sin(t178);
	t182 = t173 * t174;
	t175 = cos(t178);
	t181 = t173 * t175;
	t177 = qJD(4) + qJD(5);
	t180 = t177 * t174;
	t179 = t177 * t175;
	t168 = -t171 * t180 + t172 * t181;
	t167 = t171 * t179 + t172 * t182;
	t166 = t171 * t181 + t172 * t180;
	t165 = t171 * t182 - t172 * t179;
	t1 = [0, 0, 0, -t180, -t180; t166, t166, t166, t167, t167; -t168, -t168, -t168, t165, t165; 0, 0, 0, -t179, -t179; -t165, -t165, -t165, t168, t168; t167, t167, t167, t166, t166; 0, 0, 0, 0, 0; -t183, -t183, -t183, 0, 0; -t184, -t184, -t184, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end