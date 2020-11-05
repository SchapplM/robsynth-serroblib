% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRP3 (for one body)
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

function Tc_mdh = S5PRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:20
	% EndTime: 2020-11-04 19:58:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:20
	% EndTime: 2020-11-04 19:58:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(pkin(7));
	t45 = sin(pkin(7));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:20
	% EndTime: 2020-11-04 19:58:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t50 = cos(qJ(2));
	t49 = sin(qJ(2));
	t48 = cos(pkin(7));
	t47 = sin(pkin(7));
	t1 = [t48 * t50, -t48 * t49, t47, t48 * pkin(1) + t47 * pkin(5) + 0; t47 * t50, -t47 * t49, -t48, t47 * pkin(1) - t48 * pkin(5) + 0; t49, t50, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:20
	% EndTime: 2020-11-04 19:58:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t57 = -qJ(3) - pkin(5);
	t56 = cos(pkin(7));
	t55 = sin(pkin(7));
	t54 = qJ(2) + pkin(8);
	t53 = cos(t54);
	t52 = sin(t54);
	t51 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t56 * t53, -t56 * t52, t55, t56 * t51 - t55 * t57 + 0; t55 * t53, -t55 * t52, -t56, t55 * t51 + t56 * t57 + 0; t52, t53, 0, sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:20
	% EndTime: 2020-11-04 19:58:20
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t62 = sin(pkin(7));
	t65 = sin(qJ(4));
	t71 = t62 * t65;
	t66 = cos(qJ(4));
	t70 = t62 * t66;
	t63 = cos(pkin(7));
	t69 = t63 * t65;
	t68 = t63 * t66;
	t61 = qJ(2) + pkin(8);
	t59 = sin(t61);
	t60 = cos(t61);
	t67 = pkin(3) * t60 + pkin(6) * t59 + cos(qJ(2)) * pkin(2) + pkin(1);
	t64 = -qJ(3) - pkin(5);
	t1 = [t60 * t68 + t71, -t60 * t69 + t70, t63 * t59, -t62 * t64 + t67 * t63 + 0; t60 * t70 - t69, -t60 * t71 - t68, t62 * t59, t67 * t62 + t63 * t64 + 0; t59 * t66, -t59 * t65, -t60, t59 * pkin(3) - t60 * pkin(6) + sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:58:20
	% EndTime: 2020-11-04 19:58:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->22), mult. (40->24), div. (0->0), fcn. (53->8), ass. (0->16)
	t77 = sin(pkin(7));
	t81 = sin(qJ(4));
	t88 = t77 * t81;
	t82 = cos(qJ(4));
	t87 = t77 * t82;
	t78 = cos(pkin(7));
	t86 = t78 * t81;
	t85 = t78 * t82;
	t84 = pkin(4) * t81 + pkin(5) + qJ(3);
	t72 = t82 * pkin(4) + pkin(3);
	t76 = qJ(2) + pkin(8);
	t74 = sin(t76);
	t75 = cos(t76);
	t79 = -qJ(5) - pkin(6);
	t83 = t72 * t75 - t74 * t79 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t75 * t85 + t88, -t75 * t86 + t87, t78 * t74, t84 * t77 + t83 * t78 + 0; t75 * t87 - t86, -t75 * t88 - t85, t77 * t74, t83 * t77 - t84 * t78 + 0; t74 * t82, -t74 * t81, -t75, t74 * t72 + t75 * t79 + sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end