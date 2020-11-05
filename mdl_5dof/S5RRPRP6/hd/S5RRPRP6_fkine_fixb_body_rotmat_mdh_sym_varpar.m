% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:51
	% EndTime: 2020-11-04 20:33:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:51
	% EndTime: 2020-11-04 20:33:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t1 = [t50, -t49, 0, 0; t49, t50, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:51
	% EndTime: 2020-11-04 20:33:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t54 = cos(qJ(1));
	t53 = cos(qJ(2));
	t52 = sin(qJ(1));
	t51 = sin(qJ(2));
	t1 = [t54 * t53, -t54 * t51, t52, t54 * pkin(1) + t52 * pkin(6) + 0; t52 * t53, -t52 * t51, -t54, t52 * pkin(1) - t54 * pkin(6) + 0; t51, t53, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:51
	% EndTime: 2020-11-04 20:33:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t59 = -qJ(3) - pkin(6);
	t58 = qJ(2) + pkin(8);
	t57 = cos(t58);
	t56 = sin(t58);
	t55 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t61 * t57, -t61 * t56, t60, t61 * t55 - t59 * t60 + 0; t60 * t57, -t60 * t56, -t61, t60 * t55 + t61 * t59 + 0; t56, t57, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:51
	% EndTime: 2020-11-04 20:33:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t69 = sin(qJ(4));
	t71 = sin(qJ(1));
	t77 = t71 * t69;
	t72 = cos(qJ(4));
	t76 = t71 * t72;
	t73 = cos(qJ(1));
	t75 = t73 * t69;
	t74 = t73 * t72;
	t70 = sin(qJ(2));
	t68 = -qJ(3) - pkin(6);
	t67 = cos(pkin(8));
	t66 = sin(pkin(8));
	t65 = qJ(2) + pkin(8);
	t64 = cos(t65);
	t63 = sin(t65);
	t62 = (t67 * pkin(3) + t66 * pkin(7) + pkin(2)) * cos(qJ(2)) + (-t66 * pkin(3) + t67 * pkin(7)) * t70 + pkin(1);
	t1 = [t64 * t74 + t77, -t64 * t75 + t76, t73 * t63, t62 * t73 - t68 * t71 + 0; t64 * t76 - t75, -t64 * t77 - t74, t71 * t63, t62 * t71 + t73 * t68 + 0; t63 * t72, -t63 * t69, -t64, t70 * pkin(2) + t63 * pkin(3) - t64 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:51
	% EndTime: 2020-11-04 20:33:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->22), mult. (40->24), div. (0->0), fcn. (53->8), ass. (0->16)
	t85 = sin(qJ(4));
	t86 = sin(qJ(1));
	t94 = t86 * t85;
	t87 = cos(qJ(4));
	t93 = t86 * t87;
	t88 = cos(qJ(1));
	t92 = t88 * t85;
	t91 = t88 * t87;
	t90 = pkin(4) * t85 + pkin(6) + qJ(3);
	t78 = t87 * pkin(4) + pkin(3);
	t82 = qJ(2) + pkin(8);
	t80 = sin(t82);
	t81 = cos(t82);
	t83 = -qJ(5) - pkin(7);
	t89 = t78 * t81 - t80 * t83 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t81 * t91 + t94, -t81 * t92 + t93, t88 * t80, t90 * t86 + t89 * t88 + 0; t81 * t93 - t92, -t81 * t94 - t91, t86 * t80, t89 * t86 - t90 * t88 + 0; t80 * t87, -t80 * t85, -t81, t80 * t78 + t81 * t83 + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end