% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:21
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:16
	% EndTime: 2020-11-04 20:21:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:16
	% EndTime: 2020-11-04 20:21:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:16
	% EndTime: 2020-11-04 20:21:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t49 = cos(pkin(8));
	t48 = sin(pkin(8));
	t1 = [t51 * t49, -t51 * t48, t50, t51 * pkin(1) + t50 * qJ(2) + 0; t50 * t49, -t50 * t48, -t51, t50 * pkin(1) - t51 * qJ(2) + 0; t48, t49, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:16
	% EndTime: 2020-11-04 20:21:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t56 = pkin(6) + qJ(2);
	t55 = pkin(8) + qJ(3);
	t54 = cos(t55);
	t53 = sin(t55);
	t52 = cos(pkin(8)) * pkin(2) + pkin(1);
	t1 = [t58 * t54, -t58 * t53, t57, t58 * t52 + t56 * t57 + 0; t57 * t54, -t57 * t53, -t58, t57 * t52 - t58 * t56 + 0; t53, t54, 0, sin(pkin(8)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:16
	% EndTime: 2020-11-04 20:21:16
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t63 = sin(pkin(9));
	t66 = sin(qJ(1));
	t72 = t66 * t63;
	t64 = cos(pkin(9));
	t71 = t66 * t64;
	t67 = cos(qJ(1));
	t70 = t67 * t63;
	t69 = t67 * t64;
	t62 = pkin(8) + qJ(3);
	t60 = sin(t62);
	t61 = cos(t62);
	t68 = pkin(3) * t61 + qJ(4) * t60 + cos(pkin(8)) * pkin(2) + pkin(1);
	t65 = pkin(6) + qJ(2);
	t1 = [t61 * t69 + t72, -t61 * t70 + t71, t67 * t60, t65 * t66 + t68 * t67 + 0; t61 * t71 - t70, -t61 * t72 - t69, t66 * t60, -t67 * t65 + t68 * t66 + 0; t60 * t64, -t60 * t63, -t61, t60 * pkin(3) - t61 * qJ(4) + sin(pkin(8)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:16
	% EndTime: 2020-11-04 20:21:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t79 = pkin(9) + qJ(5);
	t75 = sin(t79);
	t84 = sin(qJ(1));
	t91 = t84 * t75;
	t77 = cos(t79);
	t90 = t84 * t77;
	t85 = cos(qJ(1));
	t89 = t85 * t75;
	t88 = t85 * t77;
	t87 = sin(pkin(9)) * pkin(4) + pkin(6) + qJ(2);
	t73 = cos(pkin(9)) * pkin(4) + pkin(3);
	t80 = pkin(8) + qJ(3);
	t76 = sin(t80);
	t78 = cos(t80);
	t82 = -pkin(7) - qJ(4);
	t86 = t73 * t78 - t76 * t82 + cos(pkin(8)) * pkin(2) + pkin(1);
	t1 = [t78 * t88 + t91, -t78 * t89 + t90, t85 * t76, t87 * t84 + t86 * t85 + 0; t78 * t90 - t89, -t78 * t91 - t88, t84 * t76, t86 * t84 - t87 * t85 + 0; t76 * t77, -t76 * t75, -t78, t76 * t73 + t78 * t82 + sin(pkin(8)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end