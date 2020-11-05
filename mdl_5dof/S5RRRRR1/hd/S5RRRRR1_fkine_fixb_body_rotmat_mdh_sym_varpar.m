% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:08
	% EndTime: 2020-11-04 20:48:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:08
	% EndTime: 2020-11-04 20:48:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [t41, -t40, 0, 0; t40, t41, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:08
	% EndTime: 2020-11-04 20:48:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (6->6), div. (0->0), fcn. (14->4), ass. (0->5)
	t45 = cos(qJ(1));
	t44 = cos(qJ(2));
	t43 = sin(qJ(1));
	t42 = sin(qJ(2));
	t1 = [t45 * t44, -t45 * t42, -t43, t45 * pkin(1) + 0; t43 * t44, -t43 * t42, t45, t43 * pkin(1) + 0; -t42, -t44, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:08
	% EndTime: 2020-11-04 20:48:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (17->11), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->7)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t49 = qJ(2) + qJ(3);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t51 * t48, -t51 * t47, -t50, t51 * t46 + 0; t50 * t48, -t50 * t47, t51, t50 * t46 + 0; -t47, -t48, 0, -sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:08
	% EndTime: 2020-11-04 20:48:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (29->14), mult. (12->10), div. (0->0), fcn. (20->8), ass. (0->8)
	t56 = qJ(2) + qJ(3);
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t55 = qJ(4) + t56;
	t54 = cos(t55);
	t53 = sin(t55);
	t52 = pkin(3) * cos(t56) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t58 * t54, -t58 * t53, -t57, t58 * t52 + 0; t57 * t54, -t57 * t53, t58, t57 * t52 + 0; -t53, -t54, 0, -pkin(3) * sin(t56) - sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:08
	% EndTime: 2020-11-04 20:48:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (53->18), mult. (34->22), div. (0->0), fcn. (47->10), ass. (0->14)
	t64 = sin(qJ(5));
	t65 = sin(qJ(1));
	t72 = t65 * t64;
	t66 = cos(qJ(5));
	t71 = t65 * t66;
	t67 = cos(qJ(1));
	t70 = t67 * t64;
	t69 = t67 * t66;
	t63 = qJ(2) + qJ(3);
	t62 = qJ(4) + t63;
	t60 = sin(t62);
	t61 = cos(t62);
	t68 = pkin(4) * t61 + pkin(6) * t60 + pkin(3) * cos(t63) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t61 * t69 - t72, -t61 * t70 - t71, t67 * t60, t68 * t67 + 0; t61 * t71 + t70, -t61 * t72 + t69, t65 * t60, t68 * t65 + 0; -t60 * t66, t60 * t64, t61, -t60 * pkin(4) + t61 * pkin(6) - pkin(3) * sin(t63) - sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end