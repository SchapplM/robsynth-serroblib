% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:59
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:44
	% EndTime: 2020-11-04 21:59:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:44
	% EndTime: 2020-11-04 21:59:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t1 = [t55, -t54, 0, 0; t54, t55, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:44
	% EndTime: 2020-11-04 21:59:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t59 = cos(qJ(1));
	t58 = cos(qJ(2));
	t57 = sin(qJ(1));
	t56 = sin(qJ(2));
	t1 = [t59 * t58, -t59 * t56, t57, t59 * pkin(1) + t57 * pkin(7) + 0; t57 * t58, -t57 * t56, -t59, t57 * pkin(1) - t59 * pkin(7) + 0; t56, t58, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:44
	% EndTime: 2020-11-04 21:59:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t64 = cos(qJ(1));
	t63 = cos(qJ(2));
	t62 = sin(qJ(1));
	t61 = sin(qJ(2));
	t60 = t63 * pkin(2) + t61 * qJ(3) + pkin(1);
	t1 = [t62, -t64 * t63, t64 * t61, t62 * pkin(7) + t60 * t64 + 0; -t64, -t62 * t63, t62 * t61, -t64 * pkin(7) + t60 * t62 + 0; 0, -t61, -t63, t61 * pkin(2) - t63 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:44
	% EndTime: 2020-11-04 21:59:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t66 = sin(pkin(9));
	t70 = sin(qJ(1));
	t77 = t70 * t66;
	t67 = cos(pkin(9));
	t76 = t70 * t67;
	t72 = cos(qJ(1));
	t75 = t72 * t66;
	t74 = t72 * t67;
	t73 = pkin(3) + pkin(7);
	t71 = cos(qJ(2));
	t69 = sin(qJ(2));
	t68 = pkin(2) + qJ(4);
	t65 = t69 * qJ(3) + t68 * t71 + pkin(1);
	t1 = [t69 * t75 + t76, t69 * t74 - t77, t72 * t71, t65 * t72 + t73 * t70 + 0; t69 * t77 - t74, t69 * t76 + t75, t70 * t71, t65 * t70 - t73 * t72 + 0; -t71 * t66, -t71 * t67, t69, -t71 * qJ(3) + t68 * t69 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:44
	% EndTime: 2020-11-04 21:59:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->20), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->15)
	t81 = sin(pkin(9));
	t85 = sin(qJ(1));
	t91 = t85 * t81;
	t82 = cos(pkin(9));
	t90 = t85 * t82;
	t87 = cos(qJ(1));
	t89 = t87 * t81;
	t88 = t87 * t82;
	t86 = cos(qJ(2));
	t84 = sin(qJ(2));
	t83 = pkin(2) + qJ(4);
	t80 = -t81 * pkin(4) + t82 * qJ(5) - qJ(3);
	t79 = pkin(4) * t82 + qJ(5) * t81 + pkin(3) + pkin(7);
	t78 = -t80 * t84 + t83 * t86 + pkin(1);
	t1 = [t84 * t89 + t90, t87 * t86, -t84 * t88 + t91, t78 * t87 + t79 * t85 + 0; t84 * t91 - t88, t85 * t86, -t84 * t90 - t89, t78 * t85 - t79 * t87 + 0; -t86 * t81, t84, t86 * t82, t80 * t86 + t83 * t84 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:44
	% EndTime: 2020-11-04 21:59:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (51->26), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->18)
	t110 = cos(qJ(6));
	t109 = sin(qJ(6));
	t100 = sin(qJ(2));
	t101 = sin(qJ(1));
	t108 = t100 * t101;
	t103 = cos(qJ(1));
	t107 = t100 * t103;
	t104 = pkin(4) + pkin(5);
	t98 = sin(pkin(9));
	t99 = cos(pkin(9));
	t106 = t99 * qJ(5) - t104 * t98 - qJ(3);
	t105 = qJ(5) * t98 + t104 * t99 + pkin(3) + pkin(7);
	t102 = cos(qJ(2));
	t97 = qJ(4) + pkin(2) - pkin(8);
	t94 = t109 * t99 - t110 * t98;
	t93 = -t109 * t98 - t110 * t99;
	t92 = -t106 * t100 + t97 * t102 + pkin(1);
	t1 = [-t101 * t93 - t94 * t107, -t101 * t94 + t93 * t107, -t103 * t102, t105 * t101 + t92 * t103 + 0; t93 * t103 - t94 * t108, t94 * t103 + t93 * t108, -t101 * t102, t92 * t101 - t105 * t103 + 0; t102 * t94, -t102 * t93, -t100, t97 * t100 + t106 * t102 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end