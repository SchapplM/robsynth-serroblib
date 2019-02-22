% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:16
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR7_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:16:22
% EndTime: 2019-02-22 11:16:22
% DurationCPUTime: 0.03s
% Computational Cost: add. (11->9), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t63 = t57 * t56;
t58 = cos(qJ(2));
t62 = t57 * t58;
t59 = cos(qJ(1));
t61 = t59 * t56;
t60 = t59 * t58;
t55 = cos(pkin(6));
t54 = sin(pkin(6));
t53 = -t55 * t63 + t60;
t52 = t55 * t62 + t61;
t51 = t55 * t61 + t62;
t50 = -t55 * t60 + t63;
t1 = [-t50, t53, 0, 0, 0, 0; t52, t51, 0, 0, 0, 0; 0, t54 * t56, 0, 0, 0, 0; t51, t52, 0, 0, 0, 0; -t53, t50, 0, 0, 0, 0; 0, -t54 * t58, 0, 0, 0, 0; -t59 * t54, 0, 0, 0, 0, 0; -t57 * t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
