% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:43
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:24
	% EndTime: 2020-11-04 20:43:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:24
	% EndTime: 2020-11-04 20:43:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t45 = cos(qJ(1));
	t44 = sin(qJ(1));
	t1 = [t45, -t44, 0, 0; t44, t45, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:24
	% EndTime: 2020-11-04 20:43:24
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
	% StartTime: 2020-11-04 20:43:24
	% EndTime: 2020-11-04 20:43:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t56 = pkin(7) + pkin(6);
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t53 = qJ(2) + qJ(3);
	t52 = cos(t53);
	t51 = sin(t53);
	t50 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t55 * t52, -t55 * t51, t54, t55 * t50 + t56 * t54 + 0; t54 * t52, -t54 * t51, -t55, t54 * t50 - t55 * t56 + 0; t51, t52, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:24
	% EndTime: 2020-11-04 20:43:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t60 = qJ(2) + qJ(3);
	t58 = sin(t60);
	t59 = cos(t60);
	t64 = pkin(3) * t59 + qJ(4) * t58 + cos(qJ(2)) * pkin(2) + pkin(1);
	t63 = pkin(7) + pkin(6);
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [t61, -t62 * t59, t62 * t58, t63 * t61 + t62 * t64 + 0; -t62, -t61 * t59, t61 * t58, t61 * t64 - t62 * t63 + 0; 0, -t58, -t59, t58 * pkin(3) - t59 * qJ(4) + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:24
	% EndTime: 2020-11-04 20:43:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->22), mult. (37->25), div. (0->0), fcn. (50->10), ass. (0->18)
	t70 = sin(qJ(5));
	t73 = sin(qJ(1));
	t81 = t73 * t70;
	t74 = cos(qJ(5));
	t80 = t73 * t74;
	t76 = cos(qJ(1));
	t79 = t76 * t70;
	t78 = t76 * t74;
	t77 = pkin(3) + pkin(8);
	t75 = cos(qJ(3));
	t72 = sin(qJ(2));
	t71 = sin(qJ(3));
	t69 = qJ(2) + qJ(3);
	t68 = pkin(4) + pkin(6) + pkin(7);
	t67 = cos(t69);
	t66 = sin(t69);
	t65 = (qJ(4) * t71 + t77 * t75 + pkin(2)) * cos(qJ(2)) + pkin(1) + (qJ(4) * t75 - t71 * t77) * t72;
	t1 = [t66 * t79 + t80, t66 * t78 - t81, t76 * t67, t65 * t76 + t68 * t73 + 0; t66 * t81 - t78, t66 * t80 + t79, t73 * t67, t65 * t73 - t68 * t76 + 0; -t67 * t70, -t67 * t74, t66, t72 * pkin(2) - t67 * qJ(4) + t77 * t66 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end