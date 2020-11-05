% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:35
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:12
	% EndTime: 2020-11-04 20:35:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:12
	% EndTime: 2020-11-04 20:35:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t1 = [t42, -t41, 0, 0; t41, t42, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:12
	% EndTime: 2020-11-04 20:35:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t46 = cos(qJ(1));
	t45 = cos(qJ(2));
	t44 = sin(qJ(1));
	t43 = sin(qJ(2));
	t1 = [t46 * t45, -t46 * t43, t44, t46 * pkin(1) + t44 * pkin(6) + 0; t44 * t45, -t44 * t43, -t46, t44 * pkin(1) - t46 * pkin(6) + 0; t43, t45, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:12
	% EndTime: 2020-11-04 20:35:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t51 = cos(qJ(1));
	t50 = cos(qJ(2));
	t49 = sin(qJ(1));
	t48 = sin(qJ(2));
	t47 = t50 * pkin(2) + t48 * qJ(3) + pkin(1);
	t1 = [t49, -t51 * t50, t51 * t48, t49 * pkin(6) + t47 * t51 + 0; -t51, -t49 * t50, t49 * t48, -t51 * pkin(6) + t47 * t49 + 0; 0, -t48, -t50, t48 * pkin(2) - t50 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:12
	% EndTime: 2020-11-04 20:35:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t53 = sin(qJ(4));
	t55 = sin(qJ(1));
	t64 = t55 * t53;
	t56 = cos(qJ(4));
	t63 = t55 * t56;
	t58 = cos(qJ(1));
	t62 = t58 * t53;
	t61 = t58 * t56;
	t60 = pkin(2) + pkin(7);
	t59 = pkin(3) + pkin(6);
	t57 = cos(qJ(2));
	t54 = sin(qJ(2));
	t52 = t54 * qJ(3) + t60 * t57 + pkin(1);
	t1 = [t54 * t62 + t63, t54 * t61 - t64, t58 * t57, t52 * t58 + t59 * t55 + 0; t54 * t64 - t61, t54 * t63 + t62, t55 * t57, t52 * t55 - t59 * t58 + 0; -t57 * t53, -t57 * t56, t54, -t57 * qJ(3) + t60 * t54 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:35:12
	% EndTime: 2020-11-04 20:35:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->20), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->15)
	t68 = sin(qJ(4));
	t70 = sin(qJ(1));
	t78 = t70 * t68;
	t71 = cos(qJ(4));
	t77 = t70 * t71;
	t73 = cos(qJ(1));
	t76 = t73 * t68;
	t75 = t73 * t71;
	t74 = pkin(2) + pkin(7);
	t72 = cos(qJ(2));
	t69 = sin(qJ(2));
	t67 = -t68 * pkin(4) + t71 * qJ(5) - qJ(3);
	t66 = pkin(4) * t71 + qJ(5) * t68 + pkin(3) + pkin(6);
	t65 = -t67 * t69 + t74 * t72 + pkin(1);
	t1 = [t69 * t76 + t77, t73 * t72, -t69 * t75 + t78, t65 * t73 + t66 * t70 + 0; t69 * t78 - t75, t70 * t72, -t69 * t77 - t76, t65 * t70 - t66 * t73 + 0; -t72 * t68, t69, t72 * t71, t67 * t72 + t74 * t69 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end