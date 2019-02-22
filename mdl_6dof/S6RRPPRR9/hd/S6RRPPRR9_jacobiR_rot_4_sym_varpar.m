% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR9
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
% Datum: 2019-02-22 11:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR9_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:17:52
% EndTime: 2019-02-22 11:17:52
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->6), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
t60 = sin(qJ(2));
t61 = sin(qJ(1));
t67 = t61 * t60;
t62 = cos(qJ(2));
t66 = t61 * t62;
t63 = cos(qJ(1));
t65 = t63 * t60;
t64 = t63 * t62;
t59 = cos(pkin(6));
t58 = sin(pkin(6));
t57 = -t59 * t67 + t64;
t56 = t59 * t66 + t65;
t55 = t59 * t65 + t66;
t54 = t59 * t64 - t67;
t1 = [t63 * t58, 0, 0, 0, 0, 0; t61 * t58, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t54, t57, 0, 0, 0, 0; t56, t55, 0, 0, 0, 0; 0, t58 * t60, 0, 0, 0, 0; -t55, -t56, 0, 0, 0, 0; t57, t54, 0, 0, 0, 0; 0, t58 * t62, 0, 0, 0, 0;];
JR_rot  = t1;
