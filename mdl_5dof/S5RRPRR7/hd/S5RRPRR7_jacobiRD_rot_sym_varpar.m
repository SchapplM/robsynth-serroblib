% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:53
	% EndTime: 2019-12-29 19:04:53
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:51
	% EndTime: 2019-12-29 19:04:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:53
	% EndTime: 2019-12-29 19:04:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0, 0; -t44, -t44, 0, 0, 0; 0, 0, 0, 0, 0; t44, t44, 0, 0, 0; -t49, -t49, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:58
	% EndTime: 2019-12-29 19:04:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t37 = qJ(1) + qJ(2);
	t36 = qJD(1) + qJD(2);
	t34 = t36 * cos(t37);
	t33 = t36 * sin(t37);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t34, t34, 0, 0, 0; t33, t33, 0, 0, 0; 0, 0, 0, 0, 0; -t33, -t33, 0, 0, 0; t34, t34, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:52
	% EndTime: 2019-12-29 19:04:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t70 = qJ(1) + qJ(2);
	t67 = sin(t70);
	t69 = qJD(1) + qJD(2);
	t78 = t69 * t67;
	t68 = cos(t70);
	t77 = t69 * t68;
	t71 = sin(qJ(4));
	t76 = t69 * t71;
	t72 = cos(qJ(4));
	t75 = t69 * t72;
	t74 = qJD(4) * t71;
	t73 = qJD(4) * t72;
	t66 = -t67 * t74 + t68 * t75;
	t65 = t67 * t73 + t68 * t76;
	t64 = t67 * t75 + t68 * t74;
	t63 = -t67 * t76 + t68 * t73;
	t1 = [t63, t63, 0, t66, 0; t65, t65, 0, t64, 0; 0, 0, 0, -t73, 0; -t64, -t64, 0, -t65, 0; t66, t66, 0, t63, 0; 0, 0, 0, t74, 0; -t77, -t77, 0, 0, 0; -t78, -t78, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:59
	% EndTime: 2019-12-29 19:04:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (134->18), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t117 = qJ(4) + qJ(5);
	t111 = sin(t117);
	t115 = qJD(4) + qJD(5);
	t110 = t115 * t111;
	t113 = cos(t117);
	t121 = t115 * t113;
	t118 = qJ(1) + qJ(2);
	t112 = sin(t118);
	t116 = qJD(1) + qJD(2);
	t120 = t116 * t112;
	t114 = cos(t118);
	t119 = t116 * t114;
	t109 = -t112 * t110 + t113 * t119;
	t108 = t111 * t119 + t112 * t121;
	t107 = t114 * t110 + t113 * t120;
	t106 = -t111 * t120 + t114 * t121;
	t1 = [t106, t106, 0, t109, t109; t108, t108, 0, t107, t107; 0, 0, 0, -t121, -t121; -t107, -t107, 0, -t108, -t108; t109, t109, 0, t106, t106; 0, 0, 0, t110, t110; -t119, -t119, 0, 0, 0; -t120, -t120, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end