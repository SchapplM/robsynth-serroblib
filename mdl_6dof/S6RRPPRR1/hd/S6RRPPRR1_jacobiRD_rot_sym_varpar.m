% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(10);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = t47 * t53 - t48 * t54;
	t45 = t47 * t54 + t48 * t53;
	t44 = t47 * t52 + t48 * t55;
	t43 = t47 * t55 - t48 * t52;
	t1 = [t46, t43, 0, 0, 0, 0; -t44, -t45, 0, 0, 0, 0; 0, -qJD(2) * t47, 0, 0, 0, 0; t45, t44, 0, 0, 0, 0; t43, t46, 0, 0, 0, 0; 0, -qJD(2) * t48, 0, 0, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:46
	% EndTime: 2019-10-10 09:35:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (28->9), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t176 = sin(qJ(1));
	t181 = qJD(1) * t176;
	t177 = cos(qJ(1));
	t180 = qJD(1) * t177;
	t179 = qJD(2) * t176;
	t178 = qJD(2) * t177;
	t175 = qJ(2) + pkin(10);
	t174 = cos(t175);
	t173 = sin(t175);
	t172 = -t173 * t179 + t174 * t180;
	t171 = -t173 * t180 - t174 * t179;
	t170 = -t173 * t178 - t174 * t181;
	t169 = t173 * t181 - t174 * t178;
	t1 = [-t172, t169, 0, 0, 0, 0; t170, t171, 0, 0, 0, 0; 0, -qJD(2) * t173, 0, 0, 0, 0; -t181, 0, 0, 0, 0, 0; t180, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t171, t170, 0, 0, 0, 0; -t169, t172, 0, 0, 0, 0; 0, qJD(2) * t174, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (161->17), mult. (250->18), div. (0->0), fcn. (250->6), ass. (0->19)
	t140 = qJD(2) - qJD(5);
	t124 = qJ(2) + pkin(10);
	t122 = sin(t124);
	t123 = cos(t124);
	t125 = sin(qJ(5));
	t127 = cos(qJ(5));
	t141 = -t122 * t127 + t123 * t125;
	t142 = t140 * t141;
	t131 = t122 * t125 + t123 * t127;
	t117 = t140 * t131;
	t130 = qJD(1) * t141;
	t129 = qJD(1) * t131;
	t128 = cos(qJ(1));
	t126 = sin(qJ(1));
	t115 = t126 * t142 + t128 * t129;
	t114 = -t126 * t117 + t128 * t130;
	t113 = t126 * t129 - t128 * t142;
	t112 = t117 * t128 + t126 * t130;
	t1 = [-t115, -t112, 0, 0, t112, 0; -t113, t114, 0, 0, -t114, 0; 0, t142, 0, 0, -t142, 0; t114, -t113, 0, 0, t113, 0; t112, t115, 0, 0, -t115, 0; 0, t117, 0, 0, -t117, 0; qJD(1) * t126, 0, 0, 0, 0, 0; -qJD(1) * t128, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:48
	% EndTime: 2019-10-10 09:35:48
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (402->47), mult. (634->69), div. (0->0), fcn. (656->8), ass. (0->43)
	t453 = qJ(2) + pkin(10);
	t451 = sin(t453);
	t452 = cos(t453);
	t455 = sin(qJ(5));
	t458 = cos(qJ(5));
	t492 = -t451 * t458 + t452 * t455;
	t484 = qJD(2) * t455;
	t495 = -t492 * qJD(5) + t452 * t484;
	t432 = t451 * t455 + t452 * t458;
	t459 = cos(qJ(1));
	t431 = t432 * t459;
	t483 = qJD(2) * t458;
	t427 = -t451 * t483 + t495;
	t494 = t432 * qJD(5) - t451 * t484;
	t493 = t492 * qJD(1);
	t456 = sin(qJ(1));
	t464 = qJD(1) * t432;
	t475 = t456 * t483;
	t425 = -t451 * t475 + t495 * t456 + t459 * t464;
	t429 = t432 * t456;
	t454 = sin(qJ(6));
	t457 = cos(qJ(6));
	t486 = qJD(1) * t456;
	t491 = (t429 * t457 + t454 * t459) * qJD(6) + t425 * t454 + t457 * t486;
	t485 = qJD(1) * t459;
	t481 = qJD(6) * t454;
	t480 = qJD(6) * t457;
	t423 = t493 * t456 + (qJD(2) - qJD(5)) * t431;
	t430 = t492 * t459;
	t470 = -t423 * t454 + t430 * t480;
	t469 = t423 * t457 + t430 * t481;
	t424 = -t452 * t475 + t494 * t456 + t493 * t459;
	t428 = t492 * t456;
	t468 = t424 * t454 + t428 * t480;
	t467 = -t424 * t457 + t428 * t481;
	t466 = t427 * t454 + t432 * t480;
	t465 = -t427 * t457 + t432 * t481;
	t460 = -t425 * t457 + t454 * t486 + (t429 * t454 - t457 * t459) * qJD(6);
	t426 = -t452 * t483 + t494;
	t422 = -t427 * t459 + t456 * t464;
	t421 = -t454 * t485 - t422 * t457 + (-t431 * t454 - t456 * t457) * qJD(6);
	t420 = -t457 * t485 + t422 * t454 + (-t431 * t457 + t454 * t456) * qJD(6);
	t1 = [t460, -t469, 0, 0, t469, t420; t421, -t467, 0, 0, t467, -t491; 0, -t465, 0, 0, t465, t426 * t454 + t480 * t492; t491, -t470, 0, 0, t470, -t421; t420, -t468, 0, 0, t468, t460; 0, -t466, 0, 0, t466, t426 * t457 - t481 * t492; -t424, t422, 0, 0, -t422, 0; -t423, -t425, 0, 0, t425, 0; 0, t426, 0, 0, -t426, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end