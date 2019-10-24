% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR2
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

function JRD_rot = S5RRRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
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
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
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
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
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
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (63->5), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t43 = qJD(1) + qJD(2) + qJD(3);
	t42 = qJ(1) + qJ(2) + qJ(3) + pkin(9);
	t40 = t43 * cos(t42);
	t39 = t43 * sin(t42);
	t1 = [0, 0, 0, 0, 0; t39, t39, t39, 0, 0; -t40, -t40, -t40, 0, 0; 0, 0, 0, 0, 0; t40, t40, t40, 0, 0; t39, t39, t39, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (182->20), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t133 = qJ(1) + qJ(2) + qJ(3) + pkin(9);
	t131 = sin(t133);
	t134 = qJD(1) + qJD(2) + qJD(3);
	t142 = t134 * t131;
	t132 = cos(t133);
	t141 = t134 * t132;
	t135 = sin(qJ(5));
	t140 = t134 * t135;
	t136 = cos(qJ(5));
	t139 = t134 * t136;
	t138 = qJD(5) * t135;
	t137 = qJD(5) * t136;
	t128 = -t131 * t138 + t132 * t139;
	t127 = t131 * t137 + t132 * t140;
	t126 = t131 * t139 + t132 * t138;
	t125 = t131 * t140 - t132 * t137;
	t1 = [0, 0, 0, 0, -t138; t126, t126, t126, 0, t127; -t128, -t128, -t128, 0, t125; 0, 0, 0, 0, -t137; -t125, -t125, -t125, 0, t128; t127, t127, t127, 0, t126; 0, 0, 0, 0, 0; -t141, -t141, -t141, 0, 0; -t142, -t142, -t142, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end