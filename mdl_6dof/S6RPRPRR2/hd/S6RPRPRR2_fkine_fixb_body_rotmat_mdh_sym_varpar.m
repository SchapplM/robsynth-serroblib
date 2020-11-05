% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:29
	% EndTime: 2020-11-04 21:39:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:29
	% EndTime: 2020-11-04 21:39:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:29
	% EndTime: 2020-11-04 21:39:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t62 = qJ(1) + pkin(10);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [t61, -t60, 0, cos(qJ(1)) * pkin(1) + 0; t60, t61, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:29
	% EndTime: 2020-11-04 21:39:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t67 = cos(qJ(3));
	t66 = sin(qJ(3));
	t65 = qJ(1) + pkin(10);
	t64 = cos(t65);
	t63 = sin(t65);
	t1 = [t64 * t67, -t64 * t66, t63, t64 * pkin(2) + t63 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t63 * t67, -t63 * t66, -t64, t63 * pkin(2) - t64 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t66, t67, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:29
	% EndTime: 2020-11-04 21:39:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t75 = -qJ(4) - pkin(7);
	t74 = qJ(1) + pkin(10);
	t73 = qJ(3) + pkin(11);
	t72 = cos(t74);
	t71 = cos(t73);
	t70 = sin(t74);
	t69 = sin(t73);
	t68 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t72 * t71, -t72 * t69, t70, t72 * t68 - t70 * t75 + cos(qJ(1)) * pkin(1) + 0; t70 * t71, -t70 * t69, -t72, t70 * t68 + t72 * t75 + sin(qJ(1)) * pkin(1) + 0; t69, t71, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:29
	% EndTime: 2020-11-04 21:39:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->23), mult. (35->24), div. (0->0), fcn. (48->10), ass. (0->15)
	t82 = qJ(1) + pkin(10);
	t78 = sin(t82);
	t84 = sin(qJ(5));
	t90 = t78 * t84;
	t85 = cos(qJ(5));
	t89 = t78 * t85;
	t80 = cos(t82);
	t88 = t80 * t84;
	t87 = t80 * t85;
	t81 = qJ(3) + pkin(11);
	t77 = sin(t81);
	t79 = cos(t81);
	t86 = pkin(4) * t79 + pkin(8) * t77 + cos(qJ(3)) * pkin(3) + pkin(2);
	t83 = -qJ(4) - pkin(7);
	t1 = [t79 * t87 + t90, -t79 * t88 + t89, t80 * t77, cos(qJ(1)) * pkin(1) - t78 * t83 + 0 + t86 * t80; t79 * t89 - t88, -t79 * t90 - t87, t78 * t77, sin(qJ(1)) * pkin(1) + t80 * t83 + 0 + t86 * t78; t77 * t85, -t77 * t84, -t79, t77 * pkin(4) - t79 * pkin(8) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:29
	% EndTime: 2020-11-04 21:39:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->27), mult. (42->26), div. (0->0), fcn. (55->12), ass. (0->18)
	t100 = qJ(1) + pkin(10);
	t94 = sin(t100);
	t101 = qJ(5) + qJ(6);
	t97 = sin(t101);
	t110 = t94 * t97;
	t98 = cos(t101);
	t109 = t94 * t98;
	t96 = cos(t100);
	t108 = t96 * t97;
	t107 = t96 * t98;
	t106 = pkin(5) * sin(qJ(5)) + pkin(7) + qJ(4);
	t104 = -pkin(9) - pkin(8);
	t91 = cos(qJ(5)) * pkin(5) + pkin(4);
	t99 = qJ(3) + pkin(11);
	t93 = sin(t99);
	t95 = cos(t99);
	t105 = -t104 * t93 + t91 * t95 + cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t95 * t107 + t110, -t95 * t108 + t109, t96 * t93, cos(qJ(1)) * pkin(1) + 0 + t106 * t94 + t105 * t96; t95 * t109 - t108, -t95 * t110 - t107, t94 * t93, sin(qJ(1)) * pkin(1) + 0 - t106 * t96 + t105 * t94; t93 * t98, -t93 * t97, -t95, t93 * t91 + t95 * t104 + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end