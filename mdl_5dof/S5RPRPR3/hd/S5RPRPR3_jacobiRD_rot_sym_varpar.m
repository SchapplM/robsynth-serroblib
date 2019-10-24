% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:42
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
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
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t13 = qJD(1) * cos(t12);
	t10 = qJD(1) * sin(t12);
	t1 = [0, 0, 0, 0, 0; t10, 0, 0, 0, 0; -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t13, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t27 = qJD(1) + qJD(3);
	t26 = qJ(1) + pkin(8) + qJ(3);
	t24 = t27 * cos(t26);
	t23 = t27 * sin(t26);
	t1 = [0, 0, 0, 0, 0; t23, 0, t23, 0, 0; -t24, 0, -t24, 0, 0; 0, 0, 0, 0, 0; t24, 0, t24, 0, 0; t23, 0, t23, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (44->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t75 = qJ(1) + pkin(8) + qJ(3);
	t73 = sin(t75);
	t76 = qJD(1) + qJD(3);
	t84 = t76 * t73;
	t74 = cos(t75);
	t83 = t76 * t74;
	t82 = t76 * sin(pkin(9));
	t81 = t76 * cos(pkin(9));
	t80 = t73 * t82;
	t79 = t74 * t81;
	t72 = t74 * t82;
	t71 = t73 * t81;
	t1 = [0, 0, 0, 0, 0; t71, 0, t71, 0, 0; -t79, 0, -t79, 0, 0; 0, 0, 0, 0, 0; -t80, 0, -t80, 0, 0; t72, 0, t72, 0, 0; 0, 0, 0, 0, 0; -t83, 0, -t83, 0, 0; -t84, 0, -t84, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (171->19), mult. (132->27), div. (0->0), fcn. (132->6), ass. (0->24)
	t176 = qJD(1) + qJD(3);
	t177 = sin(pkin(9));
	t190 = t176 * t177;
	t178 = cos(pkin(9));
	t179 = sin(qJ(5));
	t189 = t178 * t179;
	t180 = cos(qJ(5));
	t188 = t178 * t180;
	t187 = qJD(5) * t179;
	t186 = qJD(5) * t180;
	t175 = qJ(1) + pkin(8) + qJ(3);
	t174 = cos(t175);
	t185 = t174 * t190;
	t173 = sin(t175);
	t184 = t173 * t187;
	t183 = t174 * t186;
	t182 = -t173 * t180 + t174 * t189;
	t181 = t173 * t188 - t174 * t179;
	t170 = t173 * t190;
	t167 = -t178 * t184 - t183 + (t173 * t179 + t174 * t188) * t176;
	t166 = t181 * qJD(5) + t182 * t176;
	t165 = t182 * qJD(5) + t181 * t176;
	t164 = -t178 * t183 - t184 + (t173 * t189 + t174 * t180) * t176;
	t1 = [0, 0, 0, 0, -t177 * t186; t165, 0, t165, 0, t166; -t167, 0, -t167, 0, t164; 0, 0, 0, 0, t177 * t187; -t164, 0, -t164, 0, t167; t166, 0, t166, 0, t165; 0, 0, 0, 0, 0; t170, 0, t170, 0, 0; -t185, 0, -t185, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end