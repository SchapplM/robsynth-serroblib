% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:54
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:11
	% EndTime: 2020-11-04 19:54:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:11
	% EndTime: 2020-11-04 19:54:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t53 = cos(pkin(7));
	t52 = sin(pkin(7));
	t1 = [t53, -t52, 0, 0; t52, t53, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:11
	% EndTime: 2020-11-04 19:54:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t57 = cos(pkin(7));
	t56 = cos(pkin(8));
	t55 = sin(pkin(7));
	t54 = sin(pkin(8));
	t1 = [t57 * t56, -t57 * t54, t55, t57 * pkin(1) + t55 * qJ(2) + 0; t55 * t56, -t55 * t54, -t57, t55 * pkin(1) - t57 * qJ(2) + 0; t54, t56, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:11
	% EndTime: 2020-11-04 19:54:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t64 = pkin(5) + qJ(2);
	t63 = cos(pkin(7));
	t62 = sin(pkin(7));
	t61 = pkin(8) + qJ(3);
	t60 = cos(t61);
	t59 = sin(t61);
	t58 = cos(pkin(8)) * pkin(2) + pkin(1);
	t1 = [t63 * t60, -t63 * t59, t62, t63 * t58 + t64 * t62 + 0; t62 * t60, -t62 * t59, -t63, t62 * t58 - t63 * t64 + 0; t59, t60, 0, sin(pkin(8)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:11
	% EndTime: 2020-11-04 19:54:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t69 = sin(pkin(7));
	t72 = sin(qJ(4));
	t78 = t69 * t72;
	t73 = cos(qJ(4));
	t77 = t69 * t73;
	t70 = cos(pkin(7));
	t76 = t70 * t72;
	t75 = t70 * t73;
	t68 = pkin(8) + qJ(3);
	t66 = sin(t68);
	t67 = cos(t68);
	t74 = pkin(3) * t67 + pkin(6) * t66 + cos(pkin(8)) * pkin(2) + pkin(1);
	t71 = pkin(5) + qJ(2);
	t1 = [t67 * t75 + t78, -t67 * t76 + t77, t70 * t66, t71 * t69 + t74 * t70 + 0; t67 * t77 - t76, -t67 * t78 - t75, t69 * t66, t74 * t69 - t70 * t71 + 0; t66 * t73, -t66 * t72, -t67, t66 * pkin(3) - t67 * pkin(6) + sin(pkin(8)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:11
	% EndTime: 2020-11-04 19:54:11
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (52->24), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t87 = sin(pkin(7));
	t90 = sin(qJ(4));
	t96 = t87 * t90;
	t91 = cos(qJ(4));
	t95 = t87 * t91;
	t88 = cos(pkin(7));
	t94 = t88 * t90;
	t93 = t88 * t91;
	t86 = pkin(8) + qJ(3);
	t84 = sin(t86);
	t85 = cos(t86);
	t92 = pkin(3) * t85 + pkin(6) * t84 + cos(pkin(8)) * pkin(2) + pkin(1);
	t89 = pkin(5) + qJ(2);
	t82 = t85 * t93 + t96;
	t81 = t85 * t94 - t95;
	t80 = t85 * t95 - t94;
	t79 = t85 * t96 + t93;
	t1 = [t82, t88 * t84, t81, t82 * pkin(4) + t81 * qJ(5) + t89 * t87 + t88 * t92 + 0; t80, t87 * t84, t79, t80 * pkin(4) + t79 * qJ(5) + t87 * t92 - t88 * t89 + 0; t84 * t91, -t85, t84 * t90, sin(pkin(8)) * pkin(2) - t85 * pkin(6) + qJ(1) + 0 + (pkin(4) * t91 + qJ(5) * t90 + pkin(3)) * t84; 0, 0, 0, 1;];
	Tc_mdh = t1;
end