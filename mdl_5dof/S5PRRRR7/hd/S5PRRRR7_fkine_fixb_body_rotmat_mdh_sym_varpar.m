% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:09
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:04
	% EndTime: 2020-11-04 20:09:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:04
	% EndTime: 2020-11-04 20:09:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t48 = cos(pkin(9));
	t47 = sin(pkin(9));
	t1 = [t48, -t47, 0, 0; t47, t48, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:04
	% EndTime: 2020-11-04 20:09:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t52 = cos(qJ(2));
	t51 = sin(qJ(2));
	t50 = cos(pkin(9));
	t49 = sin(pkin(9));
	t1 = [t50 * t52, -t50 * t51, t49, pkin(1) * t50 + pkin(5) * t49 + 0; t49 * t52, -t49 * t51, -t50, pkin(1) * t49 - pkin(5) * t50 + 0; t51, t52, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:04
	% EndTime: 2020-11-04 20:09:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t55 = sin(qJ(3));
	t58 = cos(qJ(2));
	t61 = t55 * t58;
	t57 = cos(qJ(3));
	t60 = t57 * t58;
	t56 = sin(qJ(2));
	t59 = pkin(2) * t58 + pkin(6) * t56 + pkin(1);
	t54 = cos(pkin(9));
	t53 = sin(pkin(9));
	t1 = [t53 * t55 + t54 * t60, t53 * t57 - t54 * t61, t54 * t56, t53 * pkin(5) + t59 * t54 + 0; t53 * t60 - t54 * t55, -t53 * t61 - t54 * t57, t53 * t56, -t54 * pkin(5) + t59 * t53 + 0; t56 * t57, -t56 * t55, -t58, t56 * pkin(2) - t58 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:04
	% EndTime: 2020-11-04 20:09:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (37->24), div. (0->0), fcn. (50->8), ass. (0->14)
	t66 = sin(pkin(9));
	t70 = cos(qJ(2));
	t75 = t66 * t70;
	t67 = cos(pkin(9));
	t74 = t67 * t70;
	t73 = pkin(3) * sin(qJ(3)) + pkin(5);
	t62 = cos(qJ(3)) * pkin(3) + pkin(2);
	t69 = sin(qJ(2));
	t71 = pkin(7) + pkin(6);
	t72 = t62 * t70 + t69 * t71 + pkin(1);
	t65 = qJ(3) + qJ(4);
	t64 = cos(t65);
	t63 = sin(t65);
	t1 = [t66 * t63 + t64 * t74, -t63 * t74 + t66 * t64, t67 * t69, t73 * t66 + t72 * t67 + 0; -t67 * t63 + t64 * t75, -t63 * t75 - t67 * t64, t66 * t69, t72 * t66 - t73 * t67 + 0; t69 * t64, -t69 * t63, -t70, t69 * t62 - t70 * t71 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:04
	% EndTime: 2020-11-04 20:09:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->26), div. (0->0), fcn. (55->10), ass. (0->15)
	t82 = qJ(3) + qJ(4);
	t90 = pkin(5) + pkin(4) * sin(t82) + pkin(3) * sin(qJ(3));
	t83 = sin(pkin(9));
	t86 = cos(qJ(2));
	t89 = t83 * t86;
	t84 = cos(pkin(9));
	t88 = t84 * t86;
	t76 = pkin(4) * cos(t82) + cos(qJ(3)) * pkin(3) + pkin(2);
	t81 = -pkin(8) - pkin(7) - pkin(6);
	t85 = sin(qJ(2));
	t87 = t76 * t86 - t81 * t85 + pkin(1);
	t80 = qJ(5) + t82;
	t79 = cos(t80);
	t78 = sin(t80);
	t1 = [t83 * t78 + t79 * t88, -t78 * t88 + t83 * t79, t84 * t85, t90 * t83 + t87 * t84 + 0; -t84 * t78 + t79 * t89, -t78 * t89 - t84 * t79, t83 * t85, t87 * t83 - t90 * t84 + 0; t85 * t79, -t85 * t78, -t86, t85 * t76 + t86 * t81 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end