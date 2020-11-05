% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:34
	% EndTime: 2020-11-04 21:29:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:34
	% EndTime: 2020-11-04 21:29:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:34
	% EndTime: 2020-11-04 21:29:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t1 = [0, -t53, t52, t53 * pkin(1) + t52 * qJ(2) + 0; 0, -t52, -t53, t52 * pkin(1) - t53 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:34
	% EndTime: 2020-11-04 21:29:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t56 = pkin(1) + qJ(3);
	t55 = cos(pkin(9));
	t54 = sin(pkin(9));
	t1 = [t57 * t54, t57 * t55, t58, t57 * qJ(2) + t56 * t58 + 0; -t58 * t54, -t58 * t55, t57, -t58 * qJ(2) + t56 * t57 + 0; t55, -t54, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:34
	% EndTime: 2020-11-04 21:29:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t63 = pkin(9) + qJ(4);
	t62 = pkin(1) + pkin(7) + qJ(3);
	t61 = cos(t63);
	t60 = sin(t63);
	t59 = sin(pkin(9)) * pkin(3) + qJ(2);
	t1 = [t64 * t60, t64 * t61, t65, t59 * t64 + t62 * t65 + 0; -t65 * t60, -t65 * t61, t64, -t59 * t65 + t62 * t64 + 0; t61, -t60, 0, cos(pkin(9)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:34
	% EndTime: 2020-11-04 21:29:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t71 = sin(qJ(5));
	t72 = sin(qJ(1));
	t79 = t72 * t71;
	t73 = cos(qJ(5));
	t78 = t72 * t73;
	t74 = cos(qJ(1));
	t77 = t74 * t71;
	t76 = t74 * t73;
	t70 = pkin(9) + qJ(4);
	t67 = sin(t70);
	t68 = cos(t70);
	t75 = pkin(4) * t67 - pkin(8) * t68 + sin(pkin(9)) * pkin(3) + qJ(2);
	t69 = pkin(1) + pkin(7) + qJ(3);
	t1 = [t67 * t78 + t77, -t67 * t79 + t76, -t72 * t68, t69 * t74 + t75 * t72 + 0; -t67 * t76 + t79, t67 * t77 + t78, t74 * t68, t69 * t72 - t75 * t74 + 0; t68 * t73, -t68 * t71, t67, t68 * pkin(4) + t67 * pkin(8) + cos(pkin(9)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:34
	% EndTime: 2020-11-04 21:29:34
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (48->24), mult. (40->24), div. (0->0), fcn. (53->8), ass. (0->16)
	t87 = sin(qJ(5));
	t88 = sin(qJ(1));
	t96 = t88 * t87;
	t89 = cos(qJ(5));
	t95 = t88 * t89;
	t90 = cos(qJ(1));
	t94 = t90 * t87;
	t93 = t90 * t89;
	t92 = pkin(5) * t87 + pkin(1) + pkin(7) + qJ(3);
	t81 = t89 * pkin(5) + pkin(4);
	t85 = pkin(9) + qJ(4);
	t82 = sin(t85);
	t83 = cos(t85);
	t86 = -qJ(6) - pkin(8);
	t91 = t81 * t82 + t83 * t86 + sin(pkin(9)) * pkin(3) + qJ(2);
	t1 = [t82 * t95 + t94, -t82 * t96 + t93, -t88 * t83, t91 * t88 + t92 * t90 + 0; -t82 * t93 + t96, t82 * t94 + t95, t90 * t83, t92 * t88 - t91 * t90 + 0; t83 * t89, -t83 * t87, t82, t83 * t81 - t82 * t86 + cos(pkin(9)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end