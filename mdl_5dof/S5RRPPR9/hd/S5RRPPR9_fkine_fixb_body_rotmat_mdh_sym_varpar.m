% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:31
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:41
	% EndTime: 2020-11-04 20:31:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:41
	% EndTime: 2020-11-04 20:31:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:41
	% EndTime: 2020-11-04 20:31:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t42 = cos(qJ(1));
	t41 = cos(qJ(2));
	t40 = sin(qJ(1));
	t39 = sin(qJ(2));
	t1 = [t42 * t41, -t42 * t39, t40, t42 * pkin(1) + t40 * pkin(6) + 0; t40 * t41, -t40 * t39, -t42, t40 * pkin(1) - t42 * pkin(6) + 0; t39, t41, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:41
	% EndTime: 2020-11-04 20:31:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t47 = cos(qJ(1));
	t46 = cos(qJ(2));
	t45 = sin(qJ(1));
	t44 = sin(qJ(2));
	t43 = t46 * pkin(2) + t44 * qJ(3) + pkin(1);
	t1 = [t47 * t46, t45, t47 * t44, t45 * pkin(6) + t43 * t47 + 0; t45 * t46, -t47, t45 * t44, -t47 * pkin(6) + t43 * t45 + 0; t44, 0, -t46, t44 * pkin(2) - t46 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:41
	% EndTime: 2020-11-04 20:31:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->16), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->8)
	t54 = pkin(2) + pkin(3);
	t53 = cos(qJ(1));
	t52 = cos(qJ(2));
	t51 = sin(qJ(1));
	t50 = sin(qJ(2));
	t49 = pkin(6) - qJ(4);
	t48 = t50 * qJ(3) + t54 * t52 + pkin(1);
	t1 = [t53 * t50, -t53 * t52, -t51, t48 * t53 + t49 * t51 + 0; t51 * t50, -t51 * t52, t53, t48 * t51 - t49 * t53 + 0; -t52, -t50, 0, -t52 * qJ(3) + t54 * t50 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:41
	% EndTime: 2020-11-04 20:31:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->15)
	t59 = sin(qJ(5));
	t61 = sin(qJ(1));
	t68 = t61 * t59;
	t62 = cos(qJ(5));
	t67 = t61 * t62;
	t64 = cos(qJ(1));
	t66 = t64 * t59;
	t65 = t64 * t62;
	t63 = cos(qJ(2));
	t60 = sin(qJ(2));
	t58 = pkin(6) - qJ(4);
	t57 = qJ(3) + pkin(4);
	t56 = pkin(2) + pkin(3) + pkin(7);
	t55 = t56 * t63 + t57 * t60 + pkin(1);
	t1 = [t60 * t65 - t68, -t60 * t66 - t67, t64 * t63, t55 * t64 + t58 * t61 + 0; t60 * t67 + t66, -t60 * t68 + t65, t61 * t63, t55 * t61 - t58 * t64 + 0; -t63 * t62, t63 * t59, t60, t56 * t60 - t57 * t63 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end