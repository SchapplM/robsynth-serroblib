% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:23
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:32
	% EndTime: 2020-11-04 21:23:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:32
	% EndTime: 2020-11-04 21:23:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:32
	% EndTime: 2020-11-04 21:23:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [t54, 0, t53, t54 * pkin(1) + t53 * qJ(2) + 0; t53, 0, -t54, t53 * pkin(1) - t54 * qJ(2) + 0; 0, 1, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:32
	% EndTime: 2020-11-04 21:23:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t61 = pkin(1) + pkin(2);
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t58 = cos(pkin(9));
	t57 = sin(pkin(9));
	t56 = -t60 * t57 + t59 * t58;
	t55 = -t59 * t57 - t60 * t58;
	t1 = [-t55, t56, 0, t59 * qJ(2) + t61 * t60 + 0; t56, t55, 0, -t60 * qJ(2) + t61 * t59 + 0; 0, 0, -1, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:32
	% EndTime: 2020-11-04 21:23:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (26->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t70 = cos(pkin(9));
	t69 = cos(pkin(10));
	t68 = sin(pkin(9));
	t67 = sin(pkin(10));
	t65 = -t68 * pkin(3) + t70 * qJ(4) - qJ(2);
	t64 = t70 * pkin(3) + t68 * qJ(4) + pkin(1) + pkin(2);
	t63 = t71 * t68 + t72 * t70;
	t62 = t72 * t68 - t71 * t70;
	t1 = [t63 * t69, -t63 * t67, t62, t64 * t72 - t65 * t71 + 0; -t62 * t69, t62 * t67, t63, t64 * t71 + t65 * t72 + 0; -t67, -t69, 0, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:32
	% EndTime: 2020-11-04 21:23:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->21), mult. (33->18), div. (0->0), fcn. (47->8), ass. (0->14)
	t77 = cos(pkin(10)) * pkin(4) + pkin(3);
	t81 = sin(pkin(9));
	t82 = cos(pkin(9));
	t83 = qJ(4) + pkin(7);
	t86 = t77 * t81 - t83 * t82 + qJ(2);
	t85 = cos(qJ(1));
	t84 = sin(qJ(1));
	t80 = pkin(10) + qJ(5);
	t79 = cos(t80);
	t78 = sin(t80);
	t75 = t84 * t81 + t85 * t82;
	t74 = t85 * t81 - t84 * t82;
	t73 = t77 * t82 + t83 * t81 + pkin(1) + pkin(2);
	t1 = [t75 * t79, -t75 * t78, t74, t73 * t85 + t86 * t84 + 0; -t74 * t79, t74 * t78, t75, t73 * t84 - t86 * t85 + 0; -t78, -t79, 0, -sin(pkin(10)) * pkin(4) - qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:32
	% EndTime: 2020-11-04 21:23:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (65->28), mult. (71->34), div. (0->0), fcn. (98->10), ass. (0->19)
	t105 = sin(qJ(1));
	t94 = pkin(10) + qJ(5);
	t93 = cos(t94);
	t97 = sin(qJ(6));
	t104 = t93 * t97;
	t98 = cos(qJ(6));
	t103 = t93 * t98;
	t102 = cos(pkin(9));
	t92 = sin(t94);
	t101 = t93 * pkin(5) + t92 * pkin(8);
	t91 = cos(pkin(10)) * pkin(4) + pkin(3);
	t95 = sin(pkin(9));
	t96 = qJ(4) + pkin(7);
	t100 = t96 * t102 - t91 * t95 - qJ(2);
	t99 = cos(qJ(1));
	t89 = t99 * t102 + t105 * t95;
	t88 = -t105 * t102 + t99 * t95;
	t87 = t91 * t102 + t96 * t95 + pkin(1) + pkin(2);
	t1 = [t89 * t103 + t88 * t97, -t89 * t104 + t88 * t98, t89 * t92, -t100 * t105 + t101 * t89 + t87 * t99 + 0; -t88 * t103 + t89 * t97, t88 * t104 + t89 * t98, -t88 * t92, t100 * t99 - t101 * t88 + t87 * t105 + 0; -t92 * t98, t92 * t97, t93, -t92 * pkin(5) + t93 * pkin(8) - sin(pkin(10)) * pkin(4) - qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end