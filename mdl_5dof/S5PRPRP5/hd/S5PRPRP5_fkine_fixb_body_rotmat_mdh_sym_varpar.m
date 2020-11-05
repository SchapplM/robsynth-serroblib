% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:58
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:52
	% EndTime: 2020-11-04 19:58:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:52
	% EndTime: 2020-11-04 19:58:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t54 = cos(pkin(7));
	t53 = sin(pkin(7));
	t1 = [t54, -t53, 0, 0; t53, t54, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:52
	% EndTime: 2020-11-04 19:58:52
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t58 = cos(qJ(2));
	t57 = sin(qJ(2));
	t56 = cos(pkin(7));
	t55 = sin(pkin(7));
	t1 = [t56 * t58, -t56 * t57, t55, t56 * pkin(1) + t55 * pkin(5) + 0; t55 * t58, -t55 * t57, -t56, t55 * pkin(1) - t56 * pkin(5) + 0; t57, t58, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:52
	% EndTime: 2020-11-04 19:58:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t60 = sin(pkin(7));
	t64 = cos(qJ(2));
	t67 = t60 * t64;
	t62 = cos(pkin(7));
	t66 = t62 * t64;
	t63 = sin(qJ(2));
	t65 = pkin(2) * t64 + qJ(3) * t63 + pkin(1);
	t61 = cos(pkin(8));
	t59 = sin(pkin(8));
	t1 = [t60 * t59 + t61 * t66, -t59 * t66 + t60 * t61, t62 * t63, t60 * pkin(5) + t65 * t62 + 0; -t62 * t59 + t61 * t67, -t59 * t67 - t62 * t61, t60 * t63, -t62 * pkin(5) + t65 * t60 + 0; t63 * t61, -t63 * t59, -t64, t63 * pkin(2) - t64 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:53
	% EndTime: 2020-11-04 19:58:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (35->24), div. (0->0), fcn. (48->8), ass. (0->14)
	t73 = sin(pkin(7));
	t77 = cos(qJ(2));
	t80 = t73 * t77;
	t74 = cos(pkin(7));
	t79 = t74 * t77;
	t69 = cos(pkin(8)) * pkin(3) + pkin(2);
	t75 = qJ(3) + pkin(6);
	t76 = sin(qJ(2));
	t78 = t69 * t77 + t75 * t76 + pkin(1);
	t72 = pkin(8) + qJ(4);
	t71 = cos(t72);
	t70 = sin(t72);
	t68 = sin(pkin(8)) * pkin(3) + pkin(5);
	t1 = [t73 * t70 + t71 * t79, -t70 * t79 + t73 * t71, t74 * t76, t68 * t73 + t78 * t74 + 0; -t74 * t70 + t71 * t80, -t70 * t80 - t74 * t71, t73 * t76, -t68 * t74 + t78 * t73 + 0; t76 * t71, -t76 * t70, -t77, t76 * t69 - t77 * t75 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:53
	% EndTime: 2020-11-04 19:58:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->24), mult. (55->30), div. (0->0), fcn. (72->8), ass. (0->18)
	t90 = sin(pkin(7));
	t94 = cos(qJ(2));
	t97 = t90 * t94;
	t91 = cos(pkin(7));
	t96 = t91 * t94;
	t86 = cos(pkin(8)) * pkin(3) + pkin(2);
	t92 = qJ(3) + pkin(6);
	t93 = sin(qJ(2));
	t95 = t86 * t94 + t92 * t93 + pkin(1);
	t89 = pkin(8) + qJ(4);
	t88 = cos(t89);
	t87 = sin(t89);
	t85 = sin(pkin(8)) * pkin(3) + pkin(5);
	t84 = t90 * t87 + t88 * t96;
	t83 = t87 * t96 - t90 * t88;
	t82 = -t91 * t87 + t88 * t97;
	t81 = t87 * t97 + t91 * t88;
	t1 = [t84, t91 * t93, t83, t84 * pkin(4) + t83 * qJ(5) + t85 * t90 + t95 * t91 + 0; t82, t90 * t93, t81, t82 * pkin(4) + t81 * qJ(5) - t85 * t91 + t95 * t90 + 0; t93 * t88, -t94, t93 * t87, -t94 * t92 + qJ(1) + 0 + (pkin(4) * t88 + qJ(5) * t87 + t86) * t93; 0, 0, 0, 1;];
	Tc_mdh = t1;
end