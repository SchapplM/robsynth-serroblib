% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPPP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:35
	% EndTime: 2020-11-04 19:39:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:35
	% EndTime: 2020-11-04 19:39:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t1 = [t43, -t42, 0, 0; t42, t43, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:35
	% EndTime: 2020-11-04 19:39:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->18), div. (0->0), fcn. (36->6), ass. (0->12)
	t44 = sin(pkin(6));
	t48 = sin(qJ(1));
	t54 = t48 * t44;
	t46 = cos(pkin(6));
	t53 = t48 * t46;
	t49 = cos(qJ(1));
	t52 = t49 * t44;
	t51 = t49 * t46;
	t45 = sin(pkin(4));
	t50 = qJ(2) * t45;
	t47 = cos(pkin(4));
	t1 = [-t47 * t54 + t51, -t47 * t53 - t52, t48 * t45, t49 * pkin(1) + t48 * t50 + 0; t47 * t52 + t53, t47 * t51 - t54, -t49 * t45, t48 * pkin(1) - t49 * t50 + 0; t45 * t44, t45 * t46, t47, t47 * qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:35
	% EndTime: 2020-11-04 19:39:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (23->18), mult. (37->24), div. (0->0), fcn. (50->6), ass. (0->14)
	t57 = sin(pkin(6));
	t61 = sin(qJ(1));
	t68 = t61 * t57;
	t59 = cos(pkin(6));
	t67 = t61 * t59;
	t62 = cos(qJ(1));
	t66 = t62 * t57;
	t65 = t62 * t59;
	t64 = t57 * pkin(2) - t59 * qJ(3);
	t58 = sin(pkin(4));
	t60 = cos(pkin(4));
	t63 = t58 * qJ(2) - t60 * t64;
	t55 = t59 * pkin(2) + t57 * qJ(3) + pkin(1);
	t1 = [t61 * t58, t60 * t68 - t65, t60 * t67 + t66, t55 * t62 + t63 * t61 + 0; -t62 * t58, -t60 * t66 - t67, -t60 * t65 + t68, t55 * t61 - t63 * t62 + 0; t60, -t58 * t57, -t58 * t59, t60 * qJ(2) + t64 * t58 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:35
	% EndTime: 2020-11-04 19:39:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->19), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->16)
	t72 = sin(pkin(6));
	t77 = sin(qJ(1));
	t84 = t77 * t72;
	t74 = cos(pkin(6));
	t83 = t77 * t74;
	t78 = cos(qJ(1));
	t82 = t78 * t72;
	t81 = t78 * t74;
	t80 = pkin(2) + qJ(4);
	t70 = t74 * qJ(3) - t80 * t72;
	t73 = sin(pkin(4));
	t75 = cos(pkin(4));
	t76 = pkin(3) + qJ(2);
	t79 = t70 * t75 + t73 * t76;
	t69 = t72 * qJ(3) + t80 * t74 + pkin(1);
	t1 = [t77 * t73, t75 * t83 + t82, -t75 * t84 + t81, t69 * t78 + t79 * t77 + 0; -t78 * t73, -t75 * t81 + t84, t75 * t82 + t83, t69 * t77 - t79 * t78 + 0; t75, -t73 * t74, t73 * t72, -t70 * t73 + t76 * t75 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end