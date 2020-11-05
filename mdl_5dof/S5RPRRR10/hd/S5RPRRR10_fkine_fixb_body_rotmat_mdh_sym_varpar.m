% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:28
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:04
	% EndTime: 2020-11-04 20:28:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:04
	% EndTime: 2020-11-04 20:28:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:04
	% EndTime: 2020-11-04 20:28:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t51 = cos(pkin(9));
	t50 = sin(pkin(9));
	t1 = [t53 * t51, -t53 * t50, t52, t53 * pkin(1) + t52 * qJ(2) + 0; t52 * t51, -t52 * t50, -t53, t52 * pkin(1) - t53 * qJ(2) + 0; t50, t51, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:04
	% EndTime: 2020-11-04 20:28:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t58 = pkin(6) + qJ(2);
	t57 = pkin(9) + qJ(3);
	t56 = cos(t57);
	t55 = sin(t57);
	t54 = cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t60 * t56, -t60 * t55, t59, t60 * t54 + t58 * t59 + 0; t59 * t56, -t59 * t55, -t60, t59 * t54 - t60 * t58 + 0; t55, t56, 0, sin(pkin(9)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:04
	% EndTime: 2020-11-04 20:28:04
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t66 = sin(qJ(4));
	t67 = sin(qJ(1));
	t74 = t67 * t66;
	t68 = cos(qJ(4));
	t73 = t67 * t68;
	t69 = cos(qJ(1));
	t72 = t69 * t66;
	t71 = t69 * t68;
	t64 = pkin(9) + qJ(3);
	t62 = sin(t64);
	t63 = cos(t64);
	t70 = pkin(3) * t63 + pkin(7) * t62 + cos(pkin(9)) * pkin(2) + pkin(1);
	t65 = pkin(6) + qJ(2);
	t1 = [t63 * t71 + t74, -t63 * t72 + t73, t69 * t62, t65 * t67 + t70 * t69 + 0; t63 * t73 - t72, -t63 * t74 - t71, t67 * t62, -t69 * t65 + t70 * t67 + 0; t62 * t68, -t62 * t66, -t63, t62 * pkin(3) - t63 * pkin(7) + sin(pkin(9)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:04
	% EndTime: 2020-11-04 20:28:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t82 = qJ(4) + qJ(5);
	t79 = sin(t82);
	t85 = sin(qJ(1));
	t93 = t85 * t79;
	t80 = cos(t82);
	t92 = t85 * t80;
	t86 = cos(qJ(1));
	t91 = t86 * t79;
	t90 = t86 * t80;
	t89 = pkin(4) * sin(qJ(4)) + pkin(6) + qJ(2);
	t76 = cos(qJ(4)) * pkin(4) + pkin(3);
	t81 = pkin(9) + qJ(3);
	t77 = sin(t81);
	t78 = cos(t81);
	t87 = -pkin(8) - pkin(7);
	t88 = t76 * t78 - t77 * t87 + cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t78 * t90 + t93, -t78 * t91 + t92, t86 * t77, t89 * t85 + t88 * t86 + 0; t78 * t92 - t91, -t78 * t93 - t90, t85 * t77, t88 * t85 - t89 * t86 + 0; t77 * t80, -t77 * t79, -t78, t77 * t76 + t78 * t87 + sin(pkin(9)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end