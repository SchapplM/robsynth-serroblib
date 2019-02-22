% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:16
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:16:39
% EndTime: 2019-02-22 10:16:39
% DurationCPUTime: 0.03s
% Computational Cost: add. (37->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t71 = sin(qJ(5));
t72 = sin(qJ(1));
t78 = t72 * t71;
t73 = cos(qJ(5));
t77 = t72 * t73;
t74 = cos(qJ(1));
t76 = t74 * t71;
t75 = t74 * t73;
t70 = pkin(9) + qJ(4);
t69 = cos(t70);
t68 = sin(t70);
t67 = t68 * t75 - t78;
t66 = t68 * t76 + t77;
t65 = t68 * t77 + t76;
t64 = -t68 * t78 + t75;
t1 = [t67, 0, 0, t69 * t77, t64, 0; t65, 0, 0, -t69 * t75, t66, 0; 0, 0, 0, -t68 * t73, -t69 * t71, 0; -t66, 0, 0, -t69 * t78, -t65, 0; t64, 0, 0, t69 * t76, t67, 0; 0, 0, 0, t68 * t71, -t69 * t73, 0; -t74 * t69, 0, 0, t72 * t68, 0, 0; -t72 * t69, 0, 0, -t74 * t68, 0, 0; 0, 0, 0, t69, 0, 0;];
JR_rot  = t1;
