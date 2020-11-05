% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:55
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:27
	% EndTime: 2020-11-04 21:55:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:27
	% EndTime: 2020-11-04 21:55:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:27
	% EndTime: 2020-11-04 21:55:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t63 = qJ(1) + pkin(11);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62, -t61, 0, cos(qJ(1)) * pkin(1) + 0; t61, t62, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:27
	% EndTime: 2020-11-04 21:55:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t68 = cos(qJ(3));
	t67 = sin(qJ(3));
	t66 = qJ(1) + pkin(11);
	t65 = cos(t66);
	t64 = sin(t66);
	t1 = [t65 * t68, -t65 * t67, t64, t65 * pkin(2) + t64 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t64 * t68, -t64 * t67, -t65, t64 * pkin(2) - t65 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t67, t68, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:27
	% EndTime: 2020-11-04 21:55:27
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t72 = sin(qJ(4));
	t75 = cos(qJ(3));
	t78 = t72 * t75;
	t74 = cos(qJ(4));
	t77 = t74 * t75;
	t73 = sin(qJ(3));
	t76 = pkin(3) * t75 + pkin(8) * t73 + pkin(2);
	t71 = qJ(1) + pkin(11);
	t70 = cos(t71);
	t69 = sin(t71);
	t1 = [t69 * t72 + t70 * t77, t69 * t74 - t70 * t78, t70 * t73, cos(qJ(1)) * pkin(1) + t69 * pkin(7) + 0 + t76 * t70; t69 * t77 - t70 * t72, -t69 * t78 - t70 * t74, t69 * t73, sin(qJ(1)) * pkin(1) - t70 * pkin(7) + 0 + t76 * t69; t73 * t74, -t73 * t72, -t75, t73 * pkin(3) - t75 * pkin(8) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:27
	% EndTime: 2020-11-04 21:55:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->23), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t85 = qJ(4) + qJ(5);
	t82 = sin(t85);
	t88 = cos(qJ(3));
	t93 = t82 * t88;
	t83 = cos(t85);
	t92 = t83 * t88;
	t91 = pkin(4) * sin(qJ(4)) + pkin(7);
	t79 = cos(qJ(4)) * pkin(4) + pkin(3);
	t87 = sin(qJ(3));
	t89 = -pkin(9) - pkin(8);
	t90 = t79 * t88 - t87 * t89 + pkin(2);
	t84 = qJ(1) + pkin(11);
	t81 = cos(t84);
	t80 = sin(t84);
	t1 = [t80 * t82 + t81 * t92, t80 * t83 - t81 * t93, t81 * t87, cos(qJ(1)) * pkin(1) + 0 + t91 * t80 + t90 * t81; t80 * t92 - t81 * t82, -t80 * t93 - t81 * t83, t80 * t87, sin(qJ(1)) * pkin(1) + 0 - t91 * t81 + t90 * t80; t87 * t83, -t87 * t82, -t88, t87 * t79 + t88 * t89 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:27
	% EndTime: 2020-11-04 21:55:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (81->26), mult. (44->28), div. (0->0), fcn. (57->12), ass. (0->16)
	t103 = qJ(4) + qJ(5);
	t109 = pkin(7) + pkin(5) * sin(t103) + pkin(4) * sin(qJ(4));
	t105 = cos(qJ(3));
	t101 = qJ(1) + pkin(11);
	t98 = sin(t101);
	t108 = t105 * t98;
	t99 = cos(t101);
	t107 = t105 * t99;
	t102 = -pkin(10) - pkin(9) - pkin(8);
	t104 = sin(qJ(3));
	t94 = pkin(5) * cos(t103) + cos(qJ(4)) * pkin(4) + pkin(3);
	t106 = -t102 * t104 + t105 * t94 + pkin(2);
	t100 = qJ(6) + t103;
	t97 = cos(t100);
	t96 = sin(t100);
	t1 = [t97 * t107 + t98 * t96, -t96 * t107 + t98 * t97, t99 * t104, cos(qJ(1)) * pkin(1) + 0 + t109 * t98 + t106 * t99; t97 * t108 - t99 * t96, -t96 * t108 - t99 * t97, t98 * t104, sin(qJ(1)) * pkin(1) + 0 - t109 * t99 + t106 * t98; t104 * t97, -t104 * t96, -t105, t105 * t102 + t104 * t94 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end