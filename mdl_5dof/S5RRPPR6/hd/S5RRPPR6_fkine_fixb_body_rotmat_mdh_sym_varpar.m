% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:52
	% EndTime: 2020-11-04 20:30:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:52
	% EndTime: 2020-11-04 20:30:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:52
	% EndTime: 2020-11-04 20:30:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t55 = cos(qJ(1));
	t54 = cos(qJ(2));
	t53 = sin(qJ(1));
	t52 = sin(qJ(2));
	t1 = [t55 * t54, -t55 * t52, t53, t55 * pkin(1) + t53 * pkin(6) + 0; t53 * t54, -t53 * t52, -t55, t53 * pkin(1) - t55 * pkin(6) + 0; t52, t54, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:52
	% EndTime: 2020-11-04 20:30:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t60 = -qJ(3) - pkin(6);
	t59 = qJ(2) + pkin(8);
	t58 = cos(t59);
	t57 = sin(t59);
	t56 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t62 * t58, -t62 * t57, t61, t62 * t56 - t60 * t61 + 0; t61 * t58, -t61 * t57, -t62, t61 * t56 + t62 * t60 + 0; t57, t58, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:52
	% EndTime: 2020-11-04 20:30:52
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t67 = sin(pkin(9));
	t73 = sin(qJ(1));
	t78 = t73 * t67;
	t69 = cos(pkin(9));
	t77 = t73 * t69;
	t74 = cos(qJ(1));
	t76 = t74 * t67;
	t75 = t74 * t69;
	t72 = sin(qJ(2));
	t71 = -qJ(3) - pkin(6);
	t70 = cos(pkin(8));
	t68 = sin(pkin(8));
	t66 = qJ(2) + pkin(8);
	t65 = cos(t66);
	t64 = sin(t66);
	t63 = (pkin(3) * t70 + qJ(4) * t68 + pkin(2)) * cos(qJ(2)) + (-pkin(3) * t68 + qJ(4) * t70) * t72 + pkin(1);
	t1 = [t65 * t75 + t78, -t65 * t76 + t77, t74 * t64, t63 * t74 - t71 * t73 + 0; t65 * t77 - t76, -t65 * t78 - t75, t73 * t64, t63 * t73 + t71 * t74 + 0; t64 * t69, -t64 * t67, -t65, pkin(2) * t72 + pkin(3) * t64 - qJ(4) * t65 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:52
	% EndTime: 2020-11-04 20:30:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t85 = pkin(9) + qJ(5);
	t81 = sin(t85);
	t90 = sin(qJ(1));
	t97 = t90 * t81;
	t83 = cos(t85);
	t96 = t90 * t83;
	t91 = cos(qJ(1));
	t95 = t91 * t81;
	t94 = t91 * t83;
	t93 = sin(pkin(9)) * pkin(4) + qJ(3) + pkin(6);
	t79 = cos(pkin(9)) * pkin(4) + pkin(3);
	t86 = qJ(2) + pkin(8);
	t82 = sin(t86);
	t84 = cos(t86);
	t89 = -pkin(7) - qJ(4);
	t92 = t79 * t84 - t82 * t89 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t84 * t94 + t97, -t84 * t95 + t96, t91 * t82, t93 * t90 + t92 * t91 + 0; t84 * t96 - t95, -t84 * t97 - t94, t90 * t82, t92 * t90 - t93 * t91 + 0; t82 * t83, -t82 * t81, -t84, t82 * t79 + t84 * t89 + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end