% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:52
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:49
	% EndTime: 2020-11-04 19:52:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:49
	% EndTime: 2020-11-04 19:52:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(pkin(7));
	t48 = sin(pkin(7));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:49
	% EndTime: 2020-11-04 19:52:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(pkin(7));
	t52 = cos(pkin(8));
	t51 = sin(pkin(7));
	t50 = sin(pkin(8));
	t1 = [t53 * t52, -t53 * t50, t51, t53 * pkin(1) + t51 * qJ(2) + 0; t51 * t52, -t51 * t50, -t53, t51 * pkin(1) - t53 * qJ(2) + 0; t50, t52, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:49
	% EndTime: 2020-11-04 19:52:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t60 = pkin(5) + qJ(2);
	t59 = cos(pkin(7));
	t58 = sin(pkin(7));
	t57 = pkin(8) + qJ(3);
	t56 = cos(t57);
	t55 = sin(t57);
	t54 = cos(pkin(8)) * pkin(2) + pkin(1);
	t1 = [t59 * t56, -t59 * t55, t58, t59 * t54 + t60 * t58 + 0; t58 * t56, -t58 * t55, -t59, t58 * t54 - t59 * t60 + 0; t55, t56, 0, sin(pkin(8)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:49
	% EndTime: 2020-11-04 19:52:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t64 = pkin(8) + qJ(3);
	t62 = sin(t64);
	t63 = cos(t64);
	t68 = pkin(3) * t63 + qJ(4) * t62 + cos(pkin(8)) * pkin(2) + pkin(1);
	t67 = pkin(5) + qJ(2);
	t66 = cos(pkin(7));
	t65 = sin(pkin(7));
	t1 = [t65, -t66 * t63, t66 * t62, t67 * t65 + t68 * t66 + 0; -t66, -t65 * t63, t65 * t62, t68 * t65 - t66 * t67 + 0; 0, -t62, -t63, t62 * pkin(3) - t63 * qJ(4) + sin(pkin(8)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:49
	% EndTime: 2020-11-04 19:52:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (41->26), div. (0->0), fcn. (54->10), ass. (0->17)
	t77 = sin(pkin(7));
	t80 = sin(qJ(5));
	t89 = t77 * t80;
	t82 = cos(qJ(5));
	t88 = t77 * t82;
	t79 = cos(pkin(7));
	t87 = t79 * t80;
	t86 = t79 * t82;
	t76 = sin(pkin(8));
	t78 = cos(pkin(8));
	t84 = pkin(3) + pkin(6);
	t85 = (qJ(4) * t76 + t84 * t78) * cos(qJ(3)) + (qJ(4) * t78 - t76 * t84) * sin(qJ(3)) + t78 * pkin(2) + pkin(1);
	t75 = pkin(8) + qJ(3);
	t74 = qJ(2) + pkin(4) + pkin(5);
	t73 = cos(t75);
	t72 = sin(t75);
	t1 = [t72 * t87 + t88, t72 * t86 - t89, t79 * t73, t74 * t77 + t85 * t79 + 0; t72 * t89 - t86, t72 * t88 + t87, t77 * t73, -t74 * t79 + t85 * t77 + 0; -t73 * t80, -t73 * t82, t72, t76 * pkin(2) - t73 * qJ(4) + t84 * t72 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end