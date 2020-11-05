% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:31
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:25
	% EndTime: 2020-11-04 20:31:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:25
	% EndTime: 2020-11-04 20:31:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t45 = cos(qJ(1));
	t44 = sin(qJ(1));
	t1 = [t45, -t44, 0, 0; t44, t45, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:25
	% EndTime: 2020-11-04 20:31:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t49 = cos(qJ(1));
	t48 = cos(qJ(2));
	t47 = sin(qJ(1));
	t46 = sin(qJ(2));
	t1 = [t49 * t48, -t49 * t46, t47, t49 * pkin(1) + t47 * pkin(6) + 0; t47 * t48, -t47 * t46, -t49, t47 * pkin(1) - t49 * pkin(6) + 0; t46, t48, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:25
	% EndTime: 2020-11-04 20:31:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t54 = cos(qJ(1));
	t53 = cos(qJ(2));
	t52 = sin(qJ(1));
	t51 = sin(qJ(2));
	t50 = t53 * pkin(2) + t51 * qJ(3) + pkin(1);
	t1 = [t54 * t53, t52, t54 * t51, t52 * pkin(6) + t50 * t54 + 0; t52 * t53, -t54, t52 * t51, -t54 * pkin(6) + t50 * t52 + 0; t51, 0, -t53, t51 * pkin(2) - t53 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:25
	% EndTime: 2020-11-04 20:31:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t67 = cos(pkin(8));
	t66 = sin(pkin(8));
	t65 = pkin(2) + pkin(3);
	t64 = cos(qJ(1));
	t63 = cos(qJ(2));
	t62 = sin(qJ(1));
	t61 = sin(qJ(2));
	t60 = pkin(6) - qJ(4);
	t57 = t61 * qJ(3) + t65 * t63 + pkin(1);
	t56 = t61 * t67 - t63 * t66;
	t55 = -t61 * t66 - t63 * t67;
	t1 = [-t64 * t55, t64 * t56, -t62, t57 * t64 + t60 * t62 + 0; -t62 * t55, t62 * t56, t64, t57 * t62 - t60 * t64 + 0; t56, t55, 0, -t63 * qJ(3) + t65 * t61 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:25
	% EndTime: 2020-11-04 20:31:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->18), mult. (32->18), div. (0->0), fcn. (46->8), ass. (0->14)
	t76 = pkin(8) + qJ(5);
	t81 = sin(t76);
	t80 = cos(qJ(1));
	t79 = cos(qJ(2));
	t78 = sin(qJ(1));
	t77 = sin(qJ(2));
	t75 = qJ(4) - pkin(6) + pkin(7);
	t74 = cos(t76);
	t73 = sin(pkin(8)) * pkin(4) + qJ(3);
	t72 = cos(pkin(8)) * pkin(4) + pkin(2) + pkin(3);
	t70 = t77 * t74 - t79 * t81;
	t69 = t79 * t74 + t77 * t81;
	t68 = t72 * t79 + t73 * t77 + pkin(1);
	t1 = [t80 * t69, t80 * t70, -t78, t68 * t80 - t75 * t78 + 0; t78 * t69, t78 * t70, t80, t68 * t78 + t75 * t80 + 0; t70, -t69, 0, t72 * t77 - t73 * t79 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end